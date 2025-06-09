#!/usr/bin/env python
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

# Slicemap script
import pymses, os
from pymses.analysis import ScalarOperator, Camera
from pymses.analysis.slicing import SliceMap
from pymses.analysis.visualization.image_plot_utils import *
from optparse import OptionParser
from time import time

parser = OptionParser()
parser.usage = "%prog ramses_directory ramses_output_number map_max_resolution=512"
(opts, args) = parser.parse_args()
try:
    fileDir = args[0]
    outNumber = int(args[1])
except:
    fileDir = None
    outNumber = None
try:
    mms = int(args[2])
except:
    mms = 512
ro = pymses.RamsesOutput(fileDir, outNumber)
outNumber = ro.iout
source = ro.amr_source(["rho"])
cam = Camera(center=[0.5, 0.5, 0.5], line_of_sight_axis='z', up_vector='y', region_size=[1., 1.],
             size_unit=ro.info["unit_length"], distance=0.0, far_cut_depth=0.0, map_max_size=mms)
op = ScalarOperator(lambda dset: dset["rho"], ro.info["unit_density"])


t0 = time()
# Optional CameraOctreeDatasource creation (may be faster...)
# from pymses.sources.ramses.octree import CameraOctreeDatasource
# esize = 0.5**(ro.info["levelmin"]+1)
# cod = CameraOctreeDatasource(cam, esize, source)
# source = cod.dset
# op = MaxLevelOperator()
# Compute the SliceMap
datamap = SliceMap(source, cam, op, verbose=False)
t1 = time()
print("SliceMap time = %.1f s" % (t1 - t0), " Octree levelmax read :", cam.get_required_resolution())
# mapPylab = apply_log_scale(map)
# import matplotlib.pyplot as P
# datamap.save_plot(cmap="jet")
# P.show()

# Datamap I/O
h5fname = "./SliceMap_%s.h5" % (outNumber)
datamap.save_HDF5(h5fname)
img = datamap.save_PNG("./SliceMap_%s.png" % (outNumber))
os.remove(h5fname)
