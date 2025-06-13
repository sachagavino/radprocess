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
from .. import pymses3

class Convert:
    def __init__(self):
        self.root_directory = os.path.abspath(os.sep)

    def update_pymsesrc(self, nb_grains=1, rho=True, vel=False, Bl=False, Br=False, P=False, Xi=False, phi=False, g=False):

        # Define the directory path
        pymses_directory = os.path.join(self.root_directory, ".pymses")

        # Check if the directory exists, and create it if it doesn't
        if not os.path.exists(pymses_directory):
            os.makedirs(pymses_directory)

        # Define the file path in the pymses directory
        rc_file = os.path.join(pymses_directory, "pymsesrc")

        # Create and write to the file
        f = open(rc_file,"w")

        f.write('{')
        f.write('    "Version": 1,')
        f.write('    "Multiprocessing max. nproc": 8,')
        f.write('    "RAMSES":{')
        f.write('        "ndimensions": 3,')
        f.write('        "amr_field_descr": [')
        if (rho==True):
            f.write(f'            {"__type__": "scalar_field", "__file_type__": "hydro", "name": "rho", "ivar": 0},')
        while i < nb_grains:
            f.write(f'            {"__type__": "scalar_field", "__file_type__": "hydro", "name": "dustratio%d", "ivar": %d},' % (i+1, i+11))
            i += 1
        if (vel==True):
            f.write(f'            {"__type__": "vector_field", "__file_type__": "hydro", "name": "vel", "ivars": [1, 2, 3]},')
        if (Bl==True):
            f.write(f'            {"__type__": "vector_field", "__file_type__": "hydro", "name": "Bl", "ivars": [4, 5, 6]},')
        if (Br==True):
            f.write(f'            {"__type__": "vector_field", "__file_type__": "hydro", "name": "Br", "ivars": [7, 8, 9]},')
        if (P==True):
            f.write(f'            {"__type__": "scalar_field", "__file_type__": "hydro", "name": "P", "ivar": 10},')
        if (Xi==True):
            f.write(f'            {"__type__": "multivalued_field", "__file_type__": "hydro", "name": "Xi", "ivar_first": 11, "nb_vars": 7},')
        if (phi==True):
            f.write(f'            {"__type__": "scalar_field", "__file_type__": "grav", "name": "phi", "ivar": 0},')
        if (g==True):
            f.write(f'            {"__type__": "vector_field", "__file_type__": "grav", "name": "g", "ivars": [1, 2, 3]}')
        f.write('        ]')
        f.write('    }')
        f.write('}')
        f.close()


    def to_radmc(self, nb_grains=1):
        self.rat2 = 1

    def to_polaris(self, nb_grains=1):
        self.rat3 = 1

    