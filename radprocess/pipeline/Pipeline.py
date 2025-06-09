import sys, os, inspect
import numpy as np
import pymses3 as pym


from . Model import Model
from . Radmc3d import Radmc3d


currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 


class Pipeline(Model):

    def do_ramses2radmc(self, nb_grains=1):
        self.convert = Radmc3d.convert()


    def do_ramses2polaris(self, nb_grains=1):
        print('do it')