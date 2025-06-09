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

# inclined galaxy rotation script local octree ray tracing movie
# Local octree reconstruction use : Pay attention to 
# the ngrid_max and misc.NUMBER_OF_PROCESSES_LIMIT parameters
# as you need enough RAM and CPU ressources available !!!
from time import time

total_time = time()
from pymses import RamsesOutput

print(("import RamsesOutput=", (time() - total_time)))
import os
from pymses.analysis.visualization import Camera, FractionOperator
from pymses.analysis.visualization.image_plot_utils import *
from pymses.analysis.visualization.raytracing import OctreeRayTracer
from pymses.sources.ramses.octree import CameraOctreeDatasource, CameraOctreeDataset
from optparse import OptionParser
import numpy as N

from pymses.utils import misc

misc.NUMBER_OF_PROCESSES_LIMIT = 1  # multiprocessing
# This ngrid_max needs to be big enough to fit your data box !
ngrid_max = 2e6  # 10e6 ~ 3.8GB of RAM memory
# If you don't have enough RAM on one node, use an other script with 
# the classic ray tracer or with the MPI ray tracer that can run on many nodes
right_eye = False  # right eye movie
img_start = 0  # Use this to restart a movie computation from this image
img_stop = 10000
nbImgRot = 150
nbImgRot2 = 100
zoom_factor = 2
nbImgTranslate = 100
zoom_factor2 = 2

# small shift to avoid grid alignment problems :
center = N.array([0.5000001, 0.5000001, 0.5000001])
region_size_init = N.array([.5, 0.28125])
distance_init = 10
far_cut_depth_init = 0.4
galaxy_angle_init = -30
perspectiveAngle = 30
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
i = 0.5  # We add .5 degree to avoid the 30/60/120... degree graphic bug
# inclined the galaxy by galaxy_angle (°) :
galaxy_angle = galaxy_angle_init
angle = N.pi * i / 180
axe = [N.cos(N.pi * galaxy_angle / 180) * N.cos(angle), N.cos(N.pi * galaxy_angle / 180) * N.sin(angle),
       -N.sin(N.pi * galaxy_angle / 180)]  # = Rz(angle) * (rotation matrix Ry(galaxy_angle) * (1,0,0))
up_vector = [N.sin(N.pi * galaxy_angle / 180) * N.cos(angle), N.sin(N.pi * galaxy_angle / 180) * N.sin(angle),
             -N.cos(N.pi * galaxy_angle / 180)]  # = Rz(angle) * (rotation matrix Ry(galaxy_angle) * (0,0,1))
region_size = region_size_init
distance = distance_init
far_cut_depth = far_cut_depth_init
cam = Camera(center=center, line_of_sight_axis=axe, up_vector=up_vector, region_size=region_size,
             size_unit=ro.info["unit_length"], distance=distance, far_cut_depth=far_cut_depth, map_max_size=mms,
             perspectiveAngle=perspectiveAngle)
# Octree source creation :
source = ro.amr_source(["rho"])
# extended by max(region_size) to allow rotation :
esize = 0.5 ** (ro.info["levelmin"] + 1) + max(max(region_size), distance, far_cut_depth)
camOctSource = cam.copy()
fullOctreeDataSource = CameraOctreeDatasource(camOctSource, esize, source,
                                              ngrid_max=ngrid_max, include_split_cells=True).dset
OctreeRT = OctreeRayTracer(fullOctreeDataSource)

t0 = time()
map, levelmax_map = OctreeRT.process(op, cam, rgb=False)
t1 = time()
if right_eye:
    i = 0.5  # change the name to avoid mess with real right eye images
else:
    i = 0
print(("rt total time = %.1f s" % (t1 - t0), "mms = ", mms, "max AMR read = ", \
    cam.get_required_resolution()))
save_map_HDF5(map, cam, map_name="img%s" % (i))
# ran used to fix the colormap during the movie 
# (= use first frame colormap for each frame)
ran = save_HDF5_to_img("./img%s.h5" % (i), cmap="jet", img_path="./")
os.remove(("./img%s.h5" % (i)))
if right_eye:
    os.remove(("./img%s.png" % (i)))  # delete this first left eye view,
    # as it is computed only to get the colormap reference
    istart = 0
else:
    istart = 1  # the first image is already computed for colormap init
for i in range(istart, nbImgRot):
    if img_start <= i <= img_stop:
        # We add .5 degree to avoid the 30/60/120... degree graphic bug
        angle = N.pi * (i + .5) / 180
        axe = [N.cos(N.pi * galaxy_angle / 180) * N.cos(angle), N.cos(N.pi * galaxy_angle / 180) * N.sin(angle),
               -N.sin(N.pi * galaxy_angle / 180)]  # = Rz(angle) * (rotation matrix Ry(galaxy_angle) * (1,0,0))
        up_vector = [N.sin(N.pi * galaxy_angle / 180) * N.cos(angle), N.sin(N.pi * galaxy_angle / 180) * N.sin(angle),
                     -N.cos(N.pi * galaxy_angle / 180)]  # = Rz(angle) * (rotation matrix Ry(galaxy_angle) * (0,0,1))
        cam = Camera(center=center, line_of_sight_axis=axe, up_vector=up_vector, region_size=region_size,
                     size_unit=ro.info["unit_length"], distance=distance, far_cut_depth=far_cut_depth, map_max_size=mms,
                     perspectiveAngle=perspectiveAngle)
        if right_eye: cam = cam.get_3D_right_eye_cam()
        t0 = time()
        map, levelmax_map = OctreeRT.process(op, cam, rgb=False)
        t1 = time()

        print(("rt total time = %.1f s" % (t1 - t0), "mms = ", mms, "max AMR read = ", \
            cam.get_required_resolution()))
        save_map_HDF5(map, cam, map_name="img%s" % (i))
        save_HDF5_to_img("./img%s.h5" % (i), cmap="jet", img_path="./", ran=ran)
        os.remove(("./img%s.h5" % (i)))
istart = nbImgRot
for i in range(istart, istart + nbImgRot2):
    if img_start <= i <= img_stop:
        galaxy_angle = galaxy_angle_init * (1 - (i - istart) / (nbImgRot2 * 1.))
        axe = N.array([N.cos(N.pi * galaxy_angle / 180) * N.cos(angle), N.cos(N.pi * galaxy_angle / 180) * N.sin(angle),
                       -N.sin(N.pi * galaxy_angle / 180)])  # = Rz(angle) * (rotation matrix Ry(galaxy_angle) * (1,0,0))
        up_vector = [N.sin(N.pi * galaxy_angle / 180) * N.cos(angle), N.sin(N.pi * galaxy_angle / 180) * N.sin(angle),
                     -N.cos(N.pi * galaxy_angle / 180)]  # = Rz(angle) * (rotation matrix Ry(galaxy_angle) * (0,0,1))
        region_size = region_size_init * ((1. / zoom_factor - 1) / nbImgRot * (i - istart) + 1)
        distance = distance_init * ((1. / zoom_factor - 1) / nbImgRot * (i - istart) + 1)
        far_cut_depth = far_cut_depth_init * ((1. / zoom_factor - 1) / nbImgRot * (i - istart) + 1)
        cam = Camera(center=center, line_of_sight_axis=axe, up_vector=up_vector, region_size=region_size,
                     size_unit=ro.info["unit_length"], distance=distance, far_cut_depth=far_cut_depth, map_max_size=mms,
                     perspectiveAngle=perspectiveAngle)
        if right_eye: cam = cam.get_3D_right_eye_cam()
        t0 = time()
        map, levelmax_map = OctreeRT.process(op, cam, rgb=False)
        t1 = time()

        print(("rt total time = %.1f s" % (t1 - t0), "mms = ", mms, "max AMR read = ", \
            cam.get_required_resolution()))
        save_map_HDF5(map, cam, map_name="img%s" % (i))
        save_HDF5_to_img("./img%s.h5" % (i), cmap="jet", img_path="./", ran=ran)
        os.remove(("./img%s.h5" % (i)))
i = istart + nbImgRot2 - 1
region_size_init *= (1. / zoom_factor - 1) / nbImgRot * (i - istart) + 1
distance_init *= (1. / zoom_factor - 1) / nbImgRot * (i - istart) + 1
far_cut_depth_init *= (1. / zoom_factor - 1) / nbImgRot * (i - istart) + 1
istart = nbImgRot + nbImgRot2
for i in range(istart, istart + nbImgTranslate):
    if img_start <= i <= img_stop:
        center += axe * far_cut_depth / (2. * nbImgTranslate)
        region_size = region_size_init * ((1. / zoom_factor2 - 1) / nbImgRot * (i - istart) + 1)
        distance = distance_init * ((1. / zoom_factor2 - 1) / nbImgRot * (i - istart) + 1)
        far_cut_depth = far_cut_depth_init * ((1. / zoom_factor2 - 1) / nbImgRot * (i - istart) + 1)
        cam = Camera(center=center, line_of_sight_axis=axe, up_vector=up_vector, region_size=region_size,
                     size_unit=ro.info["unit_length"], distance=distance, far_cut_depth=far_cut_depth, map_max_size=mms,
                     perspectiveAngle=perspectiveAngle)
        if right_eye: cam = cam.get_3D_right_eye_cam()
        t0 = time()
        map, levelmax_map = OctreeRT.process(op, cam, rgb=False)
        t1 = time()

        print(("rt total time = %.1f s" % (t1 - t0), "mms = ", mms, "max AMR read = ", \
            cam.get_required_resolution()))
        save_map_HDF5(map, cam, map_name="img%s" % (i))
        save_HDF5_to_img("./img%s.h5" % (i), cmap="jet", img_path="./", ran=ran)
        os.remove(("./img%s.h5" % (i)))

print(("total film time=", (time() - total_time)))
