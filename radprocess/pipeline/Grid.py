"""
_____________________________________________________________________________________________________________
file name: Grid
@author: S. Gavino for chemistry codes.
last update: jan 2026
language: PYTHON 3.12 (constrained by pymses)
short description:  class Grid for both radmc3d and polaris. 
_____________________________________________________________________________________________________________
"""

import numpy as np

class Grid:
    def __init__(self):
        self.density = []
        self.dustdensity = []
        self.temperature = []
        self.localfield = []
        self.stars = []
        self.isrf = []
        self.dust = []
        self.accretionheating = []
        self.amr_grid = []
        self.radmc_grid = []
        self.polaris_grid = []

    def add_star(self, star):
        self.stars.append(star)

    def add_isrf(self, isrf):
        self.isrf.append(isrf)

    def add_temperature(self, temperature):
        self.temperature.append(temperature)

    def add_localfield(self, localfield):
        self.localfield.append(localfield)

    def add_density(self, density):
        self.density.append(density)

    def add_radmc_grid(self, radmc_grid):
        self.radmc_grid.append(radmc_grid)

    def add_polaris_grid(self, polaris_grid):
        self.polaris_grid.append(polaris_grid)

    def add_amr_grid(self, amr_grid):
        self.amr_grid.append(amr_grid)

    def add_dust(self, dust):
        self.dust.append(dust)

    def set_wavelength_grid(self, lmin, lmax, nlam, log=False): #microns
        if log:
            self.lam = np.logspace(np.log10(lmin), np.log10(lmax), \
                    nlam)
        else:
            self.lam = np.linspace(lmin, lmax, nlam)

    def set_mcmonowavelength_grid(self, lmin, lmax, nlam, log=False):
        if log:
            self.monolam = np.logspace(np.log10(lmin), np.log10(lmax), \
                    nlam)
        else:
            self.monolam = np.linspace(lmin, lmax, nlam)