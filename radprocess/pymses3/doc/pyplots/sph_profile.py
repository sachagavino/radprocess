# -*- coding: utf-8 -*-
#   This file is part of PyMSES.
#
#   PyMSES is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   PyMSES is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with PyMSES.  If not, see <http://www.gnu.org/licenses/>.

import numpy
import pylab
import os
from pymses import RamsesOutput
from pymses.filters import RegionFilter, PointFunctionFilter
from pymses.utils.regions import Sphere
from pymses.analysis import bin_spherical
from pymses.utils import constants as C

# Halo parameters
halo_center = [ 0.567811, 0.586055, 0.559156 ]        # in box units
halo_radius = 0.00075                                  # in box units


# RamsesOutput
datadir = os.path.join(os.path.expanduser("~"), "data/Aquarius/output")
ro = RamsesOutput(datadir, 193)


# Prepare to read the mass/epoch fields only
source = ro.particle_source(["mass", "epoch"])


# Sphere region
sph = Sphere(halo_center, halo_radius)


# Filtering particles
point_dset = RegionFilter(sph, source)
dm_filter = lambda dset: dset["epoch"] == 0.0
dm_parts = PointFunctionFilter(dm_filter, point_dset)


# Profile computation
m_weight_func = lambda dset: dset["mass"]
r_bins = numpy.linspace(0.0, halo_radius, 200)

# Mass profile
# This triggers the actual reading of the particle data files from disk.
mass_profile = bin_spherical(dm_parts, halo_center, m_weight_func, r_bins, divide_by_counts=False)

# Density profile
sph_vol = 4.0/3.0 * numpy.pi * r_bins**3
shell_vol = numpy.diff(sph_vol)
rho_profile = mass_profile / shell_vol


# Plot
# Geometrical midpoint of the bins
length = ro.info["unit_length"].express(C.kpc)
bins_centers = (r_bins[1:]+r_bins[:-1])/2. * length
dens = ro.info["unit_density"].express(C.Msun/C.kpc**3)
pylab.semilogy(bins_centers, rho_profile * dens, 'r-')
pylab.xlabel("r (kpc)")
pylab.ylabel(r"$\rho$ ($M_{\odot}.kpc^{-3}$)")
pylab.title("Spherical density profile")
