"""
_____________________________________________________________________________________________________________
file name: Convert
last update: June 2025
language: PYTHON 3.8
short description: read ramses and radmc3d files. 
_____________________________________________________________________________________________________________
"""
from __future__ import absolute_import
import os, sys, inspect
import numpy as np

from ..constants.constants import mu, autocm, amu, Ggram, kb, M_sun

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

class Convert:
    def __init__(self):
        self.test = 1

    def update_pymsesrc(self):
        self.rat = 1

    def convert2radmc(self):
        self.rat2 = 1

    def convert2polaris(self):
        self.rat3 = 1

    