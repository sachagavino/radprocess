import sys, os, inspect
import numpy as np


from . Interface import Interface
from .Convert import Convert
# from . Radmc3d import Radmc3d
# from . Polaris import Polaris


currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 


class Interface(Pipeline):

    def set_pymsesrc(self, nb_grains=1):
        self.pymsesrc = Convert(nb_grains=nb_grains)

    def do_ramses2radmc(self, nb_grains=1):
        self.to_radmc = Convert(nb_grains=nb_grains)


    def do_ramses2polaris(self, nb_grains=1):
        self.to_polaris = Convert(nb_grains=nb_grains)