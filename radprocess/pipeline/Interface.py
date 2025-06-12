import sys, os, inspect
import numpy as np


from . Pipeline import Pipeline
from . Radmc3d import Radmc3d
# from . Polaris import Polaris


currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 


class Interface(Pipeline):

    def set_pymsesrc(self, nb_grains=1, rho=True, vel=True, Br=False):
        self.pymsesrc = self.convert.update_pymsesrc(nb_grains=nb_grains, rho=rho, vel=vel, Br=Br)

    def do_ramses2radmc(self, nb_grains=1):
        self.to_radmc = self.convert.to_radmc(nb_grains=nb_grains)


    def do_ramses2polaris(self, nb_grains=1):
        self.to_polaris = self.convert.to_polaris(nb_grains=nb_grains)