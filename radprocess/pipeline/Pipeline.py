import glob, os, sys, shutil 
import numpy as np

from . Convert import Convert

from ..constants.constants import au2cm, M_sun, R_sun

class Pipeline:

    def __init__(self):
        self.convert = Convert() 


    def thermal(self, nphot=1e4, run=True, write_dens=True, \
                                           write_grid=True, \
                                           write_opac=True, \
                                           write_control=True, \
                                           write_star=True, \
                                           write_wave=True, \
                                           write_mcmono=True, \
                                           write_ext=True, \
                                           **keywords):
        """ 
        Notes:
        run MC dust radiative transfer, open the resulting dust temperature as an array and computes the surface-area weigthed temperature. If run == False, user assumes the RADMC3D output files already exist.
        -----
	    """	

        self.write_radmc3d(nphot_therm=nphot, \
                           write_dens=write_dens, \
                           write_grid=write_grid, \
                           write_opac=write_opac, \
                           write_control=write_control, \
                           write_star=write_star, \
                           write_wave=write_wave, \
                           write_mcmono=write_mcmono, \
                           write_ext=write_ext, \
                           **keywords)

        if write_dens == False or write_grid == False or write_control == False or write_opac == False or write_star == False or write_wave == False:
            print('Some RADMC3D input files will not be created. Will continue... but errors can be raised if one or more required input files are missing.\n')

        if run == True:
            self.run_thermal_radmc3d(nphot=nphot, **keywords)