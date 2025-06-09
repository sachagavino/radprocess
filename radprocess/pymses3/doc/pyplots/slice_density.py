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

import os
# import matplotlib.pyplot as P
from pymses import RamsesOutput
from pymses.analysis import Camera, ScalarOperator
from pymses.analysis.slicing import SliceMap
# from pymses.analysis.plot import Plot2D

# RamsesOutput
datadir = os.path.join(os.path.expanduser("~"), "data/Aquarius/output")
ro = RamsesOutput(datadir, 193)

# AMR data source
amr = ro.amr_source(["rho"])

# Defining a Camera object
cam = Camera(center=[0.5, 0.5, 0.5], line_of_sight_axis='z', up_vector='y', region_size=[1., 1.],
             size_unit=ro.info["unit_length"], map_max_size=256, log_sensitive=True)

# Density field access operator
rho_op = ScalarOperator(lambda dset: dset["rho"], ro.info["unit_density"])

# Slice map computation
map = SliceMap(amr, cam, rho_op, z=0.4)  # create a density slice map at z=0.4 depth position

# plot_gen = Plot2D()
fig = map.save_plot()# plot_gen)
# P.show()
