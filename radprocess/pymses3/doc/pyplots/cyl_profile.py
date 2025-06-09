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
from pymses.utils.regions import Cylinder
from pymses.analysis import bin_cylindrical, point_sampling
from pymses.utils import constants as C

# Galactic cylinder parameters
gal_center = [ 0.567811, 0.586055, 0.559156 ]        # in box units
gal_radius = 0.00024132905460547268                  # in box units
gal_thickn = 0.00010238202316595811                  # in box units
gal_normal = [ -0.172935, 0.977948, -0.117099 ]      # Norm = 1

# RamsesOutput
datadir = os.path.join(os.path.expanduser("~"), "data/Aquarius/output")
ro = RamsesOutput(datadir, 193)

# Prepare to read the density field only
source = ro.amr_source(["rho"])

# Cylinder region
cyl = Cylinder(gal_center, gal_normal, gal_radius, gal_thickn)

# AMR density field point sampling
points = cyl.random_points(1.0e6) # 1M sampling points
psp = point_sampling.PointSamplerProcessor(source)
point_dset = psp.process(points)
rho_weight_func = lambda dset: dset["rho"]
r_bins = numpy.linspace(0.0, gal_radius, 200)

# Profile computation
rho_profile = bin_cylindrical(point_dset, gal_center, gal_normal,rho_weight_func, r_bins, divide_by_counts=True)

# Plot
# Geometrical midpoint of the bins
length = ro.info["unit_length"].express(C.kpc)
bins_centers = (r_bins[1:]+r_bins[:-1])/2. * length
dens = ro.info["unit_density"].express(C.H_cc)
pylab.semilogy(bins_centers, rho_profile * dens, 'r-')
pylab.xlabel("r (kpc)")
pylab.ylabel(r"$\rho$ (H/cc)")
pylab.title("Cylindrical density profile")
