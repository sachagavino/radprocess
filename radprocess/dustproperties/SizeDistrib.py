"""
_____________________________________________________________________________________________________________
file name: SizeDistrib
@author: Sacha Gavino
last update: July 2025
language: PYTHON 3.10
short description:  Creates a discretized size distribution of dust grains.
_____________________________________________________________________________________________________________
"""
import numpy as np
from .. constants.constants import mu, amu



class SizeDistrib:
    
        def __init__(self, ndust=1, amin=5.000e-03, amax=1.000e+03, \
                           d_exp=3.5, rho_m=2.5, dtogas=0.01, cst_norm=7.41e-26, ext_eff=4):            
            self.ndust = ndust          # Number of grain size populations
            self.amin = amin            # Minimum grain size in microns
            self.amax = amax            # Maximum grain size in microns
            self.d_exp = d_exp          # Dust size distribution exponent
            self.rho_m = rho_m          # Density of the dust material in g/cm^3
            self.dtogas = dtogas        # Dust-to-gas ratio
            self.cst_norm = cst_norm    # Normalization constant



        def sizes(self):
            """ A)
            Create an array with interval and averaged size values for each grain population. 
            Returns
            -------
                Numpy array (amin, amax, average). len(array) = (3, self.ndust). Units: microns
            """
            intervals = np.logspace(np.log10(self.amin), np.log10(self.amax), self.ndust+1)
            average = np.sqrt(intervals[1:] * intervals[:-1])
            sizes_param =  np.array([intervals[:-1], intervals[1:], average])  # intervals[:-1] all but last; intervals[1:] all but first. Provides the averaged, min and max values of all ranges.

            return sizes_param


        def grainmass(self):
            """ B)
            Create an array with masses of a single grain for each grain population. 
            Returns
            -------
                Numpy array. len(array) = (nb_sizes). Units: gram
            """
            a = self.sizes()
            if self.nb_sizes == 1:
                mass = ((4.*np.pi)/3.)*self.rho_m*(a[-1]*1e-4)**3
                
            if self.nb_sizes > 1:
                mass = ((4.*np.pi)/3.)*self.rho_m*(a[-1]*1e-4)**3

            return mass

        def grainmass_single(self):
            """ C)
            Create mass from single grain (rsingle). 
            Returns
            -------
                Numpy float. Units: gram
            """
            mass_single = ((4.*np.pi)/3.)*self.rho_m*(self.rsingle*1e-4)**3

            return mass_single

        def massfraction(self):
            """ C)
            Calculate the mass fraction of each grain population relative to the total dust mass of the object following a MRN distribution. 
            Returns
            -------
                Numpy array. len(array) = (nb_sizes).
            """
            a = self.sizes()

            if self.nb_sizes == 1:
                fraction = np.array([1])

            if self.nb_sizes > 1:
                cst_norm = (3./(4.*np.pi))*((self.dtogas*mu*amu)/self.rho_m)*(4 - self.d_exp)*\
                           (1/(self.amax**(4-self.d_exp) - self.amin**(4-self.d_exp)))

                mass_density = ((4.*np.pi)/3.*(4.-self.d_exp))*self.rho_m*cst_norm*(a[1]**(4.-self.d_exp) - \
	                             a[0]**(4.-self.d_exp))

                total_density = np.sum(mass_density)
                fraction = mass_density/total_density

            return fraction
