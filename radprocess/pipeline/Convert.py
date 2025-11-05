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

from radprocess import pymses3
from radprocess import radmc3d
from radprocess import ramses

class Convert:
    #def __init__(self):

        # self.nvar = ramses.read.hydro_file_descriptor(self.path_to_ramses_model)[0]
        # self.ndust = ramses.read.hydro_file_descriptor(self.path_to_ramses_model)[2]
        # self.nvar = ramses.read.hydro_file_descriptor('/Users/sachagavino/science/projects/ecogal/enygma/fiducial_test/ramses_model/maxime_model/output_00013/')[0]
        # self.ndust = ramses.read.hydro_file_descriptor('/Users/sachagavino/science/projects/ecogal/enygma/fiducial_test/ramses_model/maxime_model/output_00013/')[2]

    def update_pymsesrc(self, 
                        ndust=1,
                        rho=True, #always true
                        dustratios=True, 
                        vel=True, 
                        bl=False, 
                        br=False, 
                        p=False, 
                        xi=False, 
                        phi=False, 
                        g=False):
        # Initialize the variable i
        i = 0
        # Define the pymsesrc path
        pymses_directory = os.path.expanduser("~/.pymses")
        
        # Check if the directory exists, and create it if it doesn't
        if not os.path.exists(pymses_directory):
            os.makedirs(pymses_directory)

        # Define the file path in the pymses directory
        rc_file = os.path.join(pymses_directory, "pymsesrc")
        print(rc_file)

        # Create and write to the file
        with open(rc_file, "w") as f:
            f.write('{\n')
            f.write('    "Version": 1,\n')
            f.write('    "Multiprocessing max. nproc": 8,\n')
            f.write('    "RAMSES":{\n')
            f.write('        "ndimensions": 3,\n')
            f.write('        "amr_field_descr": [\n')
            f.write('            {"__type__": "scalar_field", "__file_type__": "hydro", "name": "rho", "ivar": 0}')
            if (dustratios==True):
                while i < ndust:
                    f.write(',\n')
                    f.write('            {"__type__": "scalar_field", "__file_type__": "hydro", "name": "dustratio%d", "ivar": %d}' % (i+1, i+11))
                    i += 1
            if (vel==True):
                f.write(',\n')
                f.write('            {"__type__": "vector_field", "__file_type__": "hydro", "name": "vel", "ivars": [1, 2, 3]}')
            if (bl==True):
                f.write(',\n')
                f.write('            {"__type__": "vector_field", "__file_type__": "hydro", "name": "Bl", "ivars": [4, 5, 6]}')
            if (br==True):
                f.write(',\n')
                f.write('            {"__type__": "vector_field", "__file_type__": "hydro", "name": "Br", "ivars": [7, 8, 9]}')
            if (p==True):
                f.write(',\n')
                f.write('            {"__type__": "scalar_field", "__file_type__": "hydro", "name": "P", "ivar": 10}')
            if (xi==True):
                f.write(',\n')
                f.write('            {"__type__": "multivalued_field", "__file_type__": "hydro", "name": "Xi", "ivar_first": 11, "nb_vars": 7}')
            if (phi==True):
                f.write(',\n')
                f.write('            {"__type__": "scalar_field", "__file_type__": "grav", "name": "phi", "ivar": 0}')
            if (g==True):
                f.write(',\n')
                f.write('            {"__type__": "vector_field", "__file_type__": "grav", "name": "g", "ivars": [1, 2, 3]}')
            f.write('\n')
            f.write('        ]\n')
            f.write('    }\n')
            f.write('}\n')



    def to_radmc(self):
        CLR_LINE =   "                                                      \r"
        cell_counter = 0
        fields = []
        i = 0
        if rho == True:
            fields.append('rho')
        if self.ndust > 0:
            while i < self.ndust:
                fields.append("dustratio%d" % (i+1))
                i += 1
        #snap = pymses3.RamsesOutput(folder,num)
        #amr = snap.amr_source(fields)
        

    def to_polaris(self):
        self.rat3 = 1

    