"""
_____________________________________________________________________________________________________________
file name: Convert
last update: June 2025
language: > PYTHON 3.8
short description: convert from ramses to radmc3d and polaris. 
_____________________________________________________________________________________________________________
"""
import os
import importlib

import numpy as np

import pymses
from pymses.utils import constants as C
from pymses.filters import CellsToPoints

from radprocess import radmc3d
from radprocess import ramses
from radprocess.pipeline.OcTree import OcTree


class Convert:

    def ramses(self, ramses_folder, ramses_num, radmc_dir, source, nb_sizes):
        importlib.reload(pymses) # Reload pymses to clear internal caches. Again a problem with Pymses.
        CLR_LINE = " " * 50 + "\r"
        cell_counter = 0
        fields = []
        i = 0
        mp = snap.info["unit_density"].express(C.g_cc) #nu*mu in gcc
        snap = pymses.RamsesOutput(ramses_folder, ramses_num)
        snap.amr_fields()
        print("the AMR source is set to: ", source)
        amr = snap.amr_source(source)
        cell_source = CellsToPoints(amr)
        cell_source.ndim
        cells = cell_source.flatten()
        output = {}
        unit_l = snap.info["unit_length"].express(C.m) # Cell lengths
        # max. number of cells
        output["dx"] = cells.get_sizes()*unit_l
        nr_of_cells = len(output["dx"])
        # Original cell positions (from 0 to 1) converted into uint length)
        output["x"] = cells.points[:,0]*unit_l
        output["y"] = cells.points[:,1]*unit_l
        output["z"] = cells.points[:,2]*unit_l
        # level of each cell
        output["level"]=np.log2(unit_l/output["dx"])
        output["dens"]  = cells["density"] # in unit of RAMSES here
        for i in range(1,nb_sizes+1):
            output["dustratio{:d}".format(i)]  = cells["dust_ratio_{:d}".format(i)]  
        if nb_sizes > 0:
            epsilon_tot = np.zeros(output["dustratio1"].shape)
            for i in range(1,nb_sizes+1):
                epsilon_tot+=output["dustratio{:d}".format(i)]
            correction_factor = 1-epsilon_tot
            rho_gas = mp*output["dens"]*correction_factor

        amr_grid = 0
    return amr_grid #radiative friendly format for inside Pipeline.py, so it is stored in Grid, then used to write radmc3d.inp files.
        

    def ramses2polaris(self):
        self.rat3 = 1

    