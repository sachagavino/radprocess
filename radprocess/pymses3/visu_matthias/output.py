import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from functools import partial
#sys.path.append("/home/tpadiole/pymses_4.1.3/")
import quantity as Q
from pymses import RamsesOutput
from pymses.utils import constants as C
from pymses.sources.ramses.filename_utils import search_valid_outputs as search_ro
from pymses.analysis.point_sampling import PointSamplingProcessor
from pymses.filters import CellsToPoints, RegionFilter
from pymses.utils.regions import Sphere
from . import myplots


# Hypotheses on Ramses version:
# - Ramses' outputs includes gravitational potential -> "gravcompat = True"

# Comment import colormaps as cmaps

# from output import Output
# import quantity as Q
# a = Output('massive_feedback_on_mass_1000_ff_sct_25_starlum_on_lmax_11/')
# a.slice_plot(Q.Er.tot,center=[0.5,0.5,0.5],size=[0.1,0.1],
#              add_sink_info=True,log_sensitive=True,vector_field=Q.vel)
# a.column_plot(typlot=Q.vel, center=[0.5,0.5,0.5],
#               size=[0.5,0.5], add_sink_info=False)


class Output(object):
    def __init__(self, pathName='.', outputNo=-1,
                 namelistName="collapse_massive.nml",
                 verbose=False):
        self._init_pathName(pathName)
        self._init_outputNo(outputNo)
        self._init_namelistName(namelistName)
        self._init_sinkDataFrame()
        self._init_massCloud()
        self._init_ramsesOutput()
        self._verbose = verbose

        self.sphere = Sphere([0.5, 0.5, 0.5], 0.25)
        self.qtyDataFrame = pd.DataFrame()
        self.slice_plot = myplots.slice(self)
        self.column_plot = myplots.column(self)
        self.IMF_plot = myplots.IMF(self)
        self.hexbin_plot = myplots.hexbin(self)
        self.hist_plot = myplots.hist(self)
        self.MeanProfile_plot = myplots.MeanProfile(self)

    @property
    def pathName(self):
        """Returns the absolute path of the simulation folder"""
        return os.path.join(self.dirName, self.folderName)

    @property
    def outputNo(self):
        """Returns the output number"""
        return self._outputNo

    @outputNo.setter
    def outputNo(self, outputNo):
        self._init_outputNo(outputNo)
        self._init_ramsesOutput()
        self._init_sinkDataFrame()
        self._init_massCloud()

    @property
    def outputName(self):
        """Returns the name of the output folder"""
        return "output_%05d" % self.outputNo

    @property
    def outputPath(self):
        """Returns the absolute path of the output folder"""
        return os.path.join(self.pathName, self.outputName)

    # Functions about namelist
    def printNamelist(self):
        """Prints on stdout the content of the namelist file"""
        ff = open(os.path.join(self.pathName, self.namelistName))
        print(ff.read())
        ff.close()

    @property
    def info(self):
        """Returns the info dict of the ramses output"""
        return self.ramsesOutput.info

    def boxSize(self, unit=C.au):
        """
        Returns the size of the simulation box in
        unit, if None returns in box unit
        """
        boxSize = 1.
        # At this point cell_size is in Box unit
        if unit is not None:
            boxSize *= self.info["unit_length"].express(unit)
        return boxSize

    def cellSize(self, level="max", unit=C.au):
        """Returns the size of a specified level in unit,
        if None returns it in box unit"""
        if level == "max":
            level = self.info["levelmax"]
        elif level == "min":
            level = self.info["levelmin"]
        cell_size = 1. / 2**level
        # At this point cell_size is in Box unit
        if unit is not None:
            cell_size *= self.info["unit_length"].express(unit)
        return cell_size

    def rAcc(self, unit=C.au):
        """Returns the accretion of the sinks in unit, if None
        returns it in box unit"""
        if("ncell_racc" in self.info):
            ncell_racc = self.info['ncell_racc']
        else:
            ncell_racc = 4
        return ncell_racc * self.cellSize("max", unit)

    def getCenter(self, center, unit=C.au):
        """Returns the center of a sink particle (specified by its Id)
        or the sink barycenter (specified by "sink_barycenter")"""
        sinkDataFrame = self.sinkDataFrame
        info = self.info
        if isinstance(center, str):
            if center == "sink_barycenter":
                center = sinkDataFrame.as_matrix(["x", "y", "z"])
                center /= info['boxlen']
                center = np.average(center, axis=0,
                                    weights=sinkDataFrame["mass"])
        elif isinstance(center, int):
            if center in sinkDataFrame.as_matrix(["Id"]):
                center = sinkDataFrame[sinkDataFrame["Id"] == center]
                center = center.as_matrix(["x", "y", "z"])[0]
                center /= info['boxlen']
            else:
                print("Sink requested (Id=%d) doesn't exist, " % center)
                print("please change the value of center")
                sys.exit(1)
        center = np.array(center)
        # At this point, center is in Box unit
        if unit is not None:
            center *= info["unit_length"].express(unit)
        return center
        
    # Functions about sink file
    @property
    def sinkPath(self):
        """Returns the path of the sink file if it exists
        Returns None otherwise"""
        sinkpath = os.path.join(self.pathName, self.outputName,
                                'sink_%05d' % self.outputNo + '.csv')
        if os.path.isfile(sinkpath):
            return sinkpath
        else:
            return None

    @property
    def massAcc(self):
        """Returns the accreted mass by the sink particles"""
        if(self.massCloud is None):
            return None
        else:
            return self.sinkDataFrame["mass"].sum() / self.massCloud

    # Functions about qtyDataFrame
    def load(self, *fields):
        """Loads a list of quantitys in the self.qtyDataFrame"""
        print("-> Loading fields ...")
        ro = self.ramsesOutput
        df = self.qtyDataFrame
        # Recherche les champs de *fields qui
        # ne sont pas encore charges dans le dataframe
        fields_to_load = []
        for field in fields:
            if field.name not in df:
                fields_to_load.append(field)
        #
        if fields_to_load:
            amr_fields_to_read = Q.get_amrfields_to_read(*fields_to_load)
            amrSource = ro.amr_source(amr_fields_to_read,
                                      grav_compat=True,
                                      verbose=self._verbose)
            amrSource = RegionFilter(self.sphere, amrSource)
            dset = CellsToPoints(amrSource).flatten()

            for field in fields_to_load:
                df[field.name] = pd.Series(field.process(ro, dset))
            print("   Fields loaded :")
            for field in fields_to_load:
                print(field.name)
            print("")
            return fields_to_load
        else:
            print("   No field to load")
            return []

    def unload(self, *fields):
        """Unloads a list of quantitys in the self.qtyDataFrame"""
        df = self.qtyDataFrame
        print("-> Unloading fields ...")
        # Unloads field in *fields if found in the _qtyDataFrame
        fields_unloaded = []
        for field in fields:
            if field.name in df:
                df.drop(labels=field.name, axis=1, inplace=True)
                fields_unloaded.append(field.name)
        if fields_unloaded:
            print("   Fields unloaded :")
            for field in fields_unloaded:
                print(field)
            print("")
        else:
            print("   No field to unload")
        if df.empty:
            self.qtyDataFrame = pd.DataFrame()

    #############################################################
    # Initialisation functions
    #############################################################

    def _init_massCloud(self):
        self.massCloud = None
        ff = open(os.path.join(self.pathName, self.namelistName))
        for line in ff.readlines():
            if('mass_c' in line):
                mass = line.split('=')[1]
                mass = mass.split('!')[0]
                self.massCloud = float(mass)
        ff.close()

    def _init_sinkDataFrame(self):
        self.sinkDataFrame = None
        print("-> Searching for sink file ...")
        if self.sinkPath is not None:
            names = ['Id', 'mass', 'dmf', 'x', 'y', 'z', 'vx', 'vy', 'vz',
                     'rot_period', 'lx', 'ly', 'lz', 'acc_rate',
                     'acc_lum', 'age', 'int_lum', 'Teff']
            self.sinkDataFrame = pd.read_csv(self.sinkPath, header=None,
                                              sep=',', names=names)
            print("OK, sink_%05d" % self.outputNo + ".csv file found")
        else:
            print("no sink file found in the folder : ")
            print(os.path.join(self.folderName, self.outputName))

    def _init_ramsesOutput(self):
        print("->")
        ro = RamsesOutput(self.pathName, self.outputNo)
        self.ramsesOutput = ro
        ro.define_amr_scalar_field("hydro", "rho", 0)
        ro.define_amr_vector_field("hydro", "vel", [1, 2, 3])
        ro.define_amr_vector_field("hydro", "Bl", [4, 5, 6])
        ro.define_amr_vector_field("hydro", "Br", [7, 8, 9])
        ro.define_amr_scalar_field("hydro", "P", 10)
        ro.define_amr_multivalued_field("hydro", "Er", 11, ro.info["ngrp"])
        ro.define_amr_vector_field("hydro", "J", [11+ro.info["ngrp"],
                                                  12+ro.info["ngrp"],
                                                  13+ro.info["ngrp"]])
        ro.define_amr_scalar_field("hydro", "eint", 14+ro.info["ngrp"])
        if ro.info["eos"]:
            ro.define_amr_scalar_field("hydro", "T", 15+ro.info["ngrp"])
        ro.define_amr_scalar_field("grav", "phi", 0)
        ro.define_amr_vector_field("grav", "g", [1, 2, 3])

    def _init_namelistName(self, namelistName):
        self.namelistName = None
        print("-> Initialising namelist ...")
        if namelistName is not None:
            if os.path.isfile(os.path.join(self.pathName, namelistName)):
                self.namelistName = namelistName
                print("OK, " + self.namelistName + " file found")
            else:
                print("ERROR, " + namelistName + " file not found ...")

        if self.namelistName is None:
            print("   Searching for namelist file in %s ..." % self.pathName)
            namelistName_list = []
            for f in os.listdir(self.pathName):
                if f.endswith(".nml"):
                    namelistName_list.append(f)
            if len(namelistName_list) != 1:
                print("ERROR, ")
                print("several namelists have been found, ")
                print("please specify the correct namelist ")
                print("at the creation of the Output object :")
                for namelistName in namelistName_list:
                    print(namelistName)
                sys.exit(1)
            else:
                self.namelistName = namelistName_list[0]
                print("OK, " + self.namelistName + " file found")

    def _init_outputNo(self, outputNo):
        print("-> Initialising output ...")
        if search_ro(self.pathName):
            if outputNo == 0:
                print("\n   Searching for the first output ...")
                self._outputNo = search_ro(self.pathName)[0]
            elif outputNo == -1:
                print("\n   Searching for the last output ...")
                self._outputNo = search_ro(self.pathName)[-1]
            elif outputNo in search_ro(self.pathName):
                self._outputNo = outputNo
            else:
                print("ERROR, output_%05d folder not found" % outputNo)
                sys.exit(1)
            print("OK, %s folder found" % self.outputName)
        else:
            print("ERROR, no valid output folder found in %s" % self.pathName)
            sys.exit(1)

    def _init_pathName(self, pathName):
        print("-> Initialising path ...")
        pathName = os.path.abspath(pathName)
        if os.path.isdir(pathName):
            self.dirName = os.path.dirname(pathName)
            self.folderName = os.path.basename(pathName)
            print("OK, %s folder found" % self.pathName)
        else:
            print("ERROR, %s folder not found" % pathName)
            sys.exit(1)
