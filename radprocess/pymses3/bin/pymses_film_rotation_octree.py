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

# basic rotation script with local octree ray tracing movie
# Local octree reconstruction use : Pay attention to 
# the ngrid_max and misc.NUMBER_OF_PROCESSES_LIMIT parameters
# as you need enough RAM and CPU ressources available !!!
from time import time

from pymses.sources.ramses.octree import CameraOctreeDatasource
from pymses.analysis import raytracing

total_time = time()
from pymses import RamsesOutput

print(("import RamsesOutput=", (time() - total_time)))
import os
from pymses.analysis import Camera, FractionOperator
from pymses.analysis.visualization.image_plot_utils import *
from optparse import OptionParser
import numpy as N

from pymses.utils import misc

misc.NUMBER_OF_PROCESSES_LIMIT = 1  # multiprocessing
# This ngrid_max needs to be big enough to fit your data box !
ngrid_max = 2e6  # 10e6 ~ 3.8GB of RAM memory
# If you don't have enough RAM on one node, use an other script with 
# the classic ray tracer or with the MPI ray tracer that can run on many nodes
right_eye = False  # to get the right eye movie here : just switch image 0 to image 1 !
img_start = 0  # Use this to restart a movie computation from this image
img_stop = 10000
nbImgRot = 18

# small shift to avoid grid alignment problems
center_init = N.array([0.500001, 0.5000001, 0.5000001])
region_size_init = N.array([.25, .140625])
distance_init = 0.2
far_cut_depth_init = 0.2
zoom_focus = N.array([0.56704, 0.58636, 0.55961])
zoom_factor = 4

parser = OptionParser()
parser.usage = "%prog ramses_directory ramses_output_number map_max_resolution=1920"
(opts, args) = parser.parse_args()
try:
    fileDir = args[0]
    outNumber = int(args[1])
    mms = int(args[2])
except:
    # "None" leads to an automatic look for
    # a RAMSES output in the current directory
    fileDir = None
    outNumber = None
    mms = 192  # 0 # 1920 <-> full HD movie
ro = RamsesOutput(fileDir, outNumber)
outNumber = ro.iout

op = FractionOperator(lambda dset: (dset["rho"] ** 2), lambda dset: (dset["rho"]))
# Colormap initialisation (ran)
# the first iteration of the following "for" loop is done separately 
# to save the colormap used for the whole movie
i = 0
angle = N.pi * i / nbImgRot + N.pi * .5 / 180  # add .5 to avoid grid alignment artifact problem
axe = [N.sin(angle), 0, N.cos(angle)]
center = center_init + (zoom_focus - center_init) * (i * 1. / nbImgRot)
region_size = region_size_init * ((1. / zoom_factor - 1) / nbImgRot * i + 1)
distance = distance_init * ((1. / zoom_factor - 1) / nbImgRot * i + 1)
far_cut_depth = far_cut_depth_init * ((1. / zoom_factor - 1) / nbImgRot * i + 1)
cam = Camera(center=center, line_of_sight_axis=axe, up_vector="y", region_size=region_size,
             size_unit=ro.info["unit_length"], distance=distance, far_cut_depth=far_cut_depth, map_max_size=mms)

# Octree source creation :
source = ro.amr_source(["rho"])
# We need to add an extension to the octree box loaded in memory,
# or not (if min(distance,far_cut_depth) > max(region_size))
# as the rotation is around the "y" axis,
# to allow rotation without empty starting point ray bugs
esize = 0.5 ** (ro.info["levelmin"] + 1) + max(region_size)
camOctSource = cam.copy()
fullOctreeDataSource = CameraOctreeDatasource(camOctSource, esize, source,
                                              ngrid_max=ngrid_max, include_split_cells=True).dset
OctreeRT = raytracing.OctreeRayTracer(fullOctreeDataSource)

t0 = time()
map, levelmax_map = OctreeRT.process(op, cam, rgb=False)
t1 = time()
print(("rt total time = %.1f s" % (t1 - t0), "mms = ", mms, "max AMR read = ", \
    cam.get_required_resolution()))
save_map_HDF5(map, cam, map_name="img%s" % (i))
# ran used to fix the colormap during the movie 
# (= use first frame colormap for each frame)
ran = save_HDF5_to_img("./img%s.h5" % (i), cmap="jet", img_path="./", ramses_output=ro)
for i in range(1, nbImgRot):
    if i >= img_start and i <= img_stop:
        # We add .5 degree to avoid a possible 30/60/120... degree graphic bug
        angle = N.pi * i / nbImgRot + N.pi * .5 / 180
        axe = [N.sin(angle), 0, N.cos(angle)]
        center = center_init + (zoom_focus - center_init) * (i * 1. / nbImgRot)
        region_size = region_size_init * ((1. / zoom_factor - 1) / nbImgRot * i + 1)
        distance = distance_init * ((1. / zoom_factor - 1) / nbImgRot * i + 1)
        far_cut_depth = far_cut_depth_init * ((1. / zoom_factor - 1) / nbImgRot * i + 1)
        cam = Camera(center=center, line_of_sight_axis=axe, up_vector="y", region_size=region_size,
                     size_unit=ro.info["unit_length"], distance=distance, far_cut_depth=far_cut_depth, map_max_size=mms)
        t0 = time()
        map, levelmax_map = OctreeRT.process(op, cam, rgb=False)
        t1 = time()
        print(("rt total time = %.1f s" % (t1 - t0), "mms = ", mms, "max AMR read = ", \
            cam.get_required_resolution()))
        save_map_HDF5(map, cam, map_name="img%s" % (i))
        save_HDF5_to_img("./img%s.h5" % (i), cmap="jet", img_path="./", ran=ran, ramses_output=ro)
        os.remove(("./img%s.h5" % (i)))
print(("total film time=", (time() - total_time)))
