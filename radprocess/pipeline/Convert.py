"""
_____________________________________________________________________________________________________________
file name: Convert
last update: June 2025
language: > PYTHON 3.8
short description: convert from ramses to radmc3d and polaris. 
_____________________________________________________________________________________________________________
"""
import os

import numpy as np

#import radprocess.pymses3 as pymses

#from radprocess.pymses.utils import constants as C
#from radprocess.pymses3.filters import CellsToPoints
from radprocess import radmc3d
from radprocess import ramses
from radprocess.pipeline.OcTree import OcTree

class Convert:
    #def __init__(self):
        # self.nvar = ramses.read.hydro_file_descriptor(self.path_to_ramses_model)[0]
        # self.ndust = ramses.read.hydro_file_descriptor(self.path_to_ramses_model)[2]
        # self.nvar = ramses.read.hydro_file_descriptor('/Users/sachagavino/science/projects/ecogal/enygma/fiducial_test/ramses_model/maxime_model/output_00013/')[0]
        # self.ndust = ramses.read.hydro_file_descriptor('/Users/sachagavino/science/projects/ecogal/enygma/fiducial_test/ramses_model/maxime_model/output_00013/')[2]


    def to_radmc(self, ramses_folder, ramses_num, radmc_dir):
        CLR_LINE =   "                                                      \r"
        cell_counter = 0
        fields = []
        i = 0
        #nvar, variables, nb_dust = ramses.read.hydro_file_descriptor(ramses_folder)
        #snap = pymses.RamsesOutput(ramses_folder, ramses_num)
        # if rho == True:
        #     fields.append('rho')
        # if self.ndust > 0:
        #     while i < self.ndust:
        #         fields.append("dustratio%d" % (i+1))
        #         i += 1
        #snap = pymses.RamsesOutput(folder,num)
        #amr = snap.amr_source(fields)
        

    def to_polaris(self):
        self.rat3 = 1

    