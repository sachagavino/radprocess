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

    def update_pymsesrc(self, nb_grains=1, rho=True, vel=True, Br=False):

        # Define the directory path
        pymses_directory = os.path.join(self.root_directory, ".pymses")

        # Check if the directory exists, and create it if it doesn't
        if not os.path.exists(pymses_directory):
            os.makedirs(pymses_directory)

        # Define the file path in the pymses directory
        file_path = os.path.join(pymses_directory, "pymsesrc")

        # Create and write to the file
        with open(file_path, "w") as file:
            file.write("This is a file in the pymses directory.")


    def to_radmc(self, nb_grains=1):
        self.rat2 = 1

    def to_polaris(self, nb_grains=1):
        self.rat3 = 1

    