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

# transfer function ray tracing map script
import pymses, os
from pymses.analysis import ScalarOperator, FractionOperator, Camera
from pymses.analysis.visualization.image_plot_utils import *
from time import time
from pymses.analysis.visualization import ColorLinesTransferFunction
from pymses.analysis.visualization.raytracing import OctreeRayTracer, RayTracer
from pymses.sources.ramses import CameraOctreeDatasource, CameraOctreeDataset
from optparse import OptionParser
from numpy import log10

t0 = time()
parser = OptionParser()
parser.usage = "%prog ramses_directory ramses_output_number map_max_resolution=512"
(opts, args) = parser.parse_args()
try:
    fileDir = args[0]
    outNumber = int(args[1])
    mms = int(args[2])
except:
    # "None" leads to an automatic look for
    # RAMSES output in the current directory
    fileDir = None
    outNumber = None
    mms = 512
ro = pymses.RamsesOutput(fileDir, outNumber)
outNumber = ro.iout
cam = Camera(center=[0.5, 0.5, 0.5], line_of_sight_axis="z", up_vector="y", region_size=[5.0E-1, 5.0E-1],
             size_unit=ro.info["unit_length"], distance=2.5E-1, far_cut_depth=2.5E-1, map_max_size=mms)
# import tables
# file= tables.openFile("camera.h5", "r")
# cam = Camera.from_HDF5(file)
# file.close()
source = ro.amr_source(["rho"])
esize = 0.5 ** (ro.info["levelmin"] + 1)

# ####################################################################################
# ####### 1) OctreeRayTracer with ColorLinesTransferFunction definition ##############
# ####################################################################################
cltf = ColorLinesTransferFunction((-5.0, 2.0))
cltf.add_line(-2.0, 0.1)
cltf.add_line(.0, 0.1)
cltf.add_line(2., 0.1)
cam.set_color_transfer_function(cltf)
# We add 1e-8 to avoid NaN and -Inf log result problems with approximative null values
op = ScalarOperator(lambda dset: log10(dset["rho"] + 1e-8), C.none)
# ##### Option A : Create OctreeRayTracer = implicit OctreeDataSource creation #######
# OctreeRT = OctreeRayTracer(ro, ["rho"])
# ##### Option B : Explicit OctreeDataSource creation (can be explicitly reused) #####
# ######## fullOctreeDataSource = Build a local octree for this aera #################
fullOctreeDataSource = CameraOctreeDatasource(cam, esize, source).dset
# Here is how to save and load local octree option for faster reuse if wanted :
# fullOctreeDataSource.save_HDF5("myDset.h5")
# fullOctreeDataSource = CameraOctreeDataset.from_HDF5("myDset.h5")
OctreeRT = OctreeRayTracer(fullOctreeDataSource)
# ########            Start OctreeRayTracer process :                 ################
img = OctreeRT.process(op, cam)
# img.show()
img.save("rt_tf_%s.png" % outNumber)
print(("rt_tf_%s.png saved" % outNumber))

# ####################################################################################
# ####### 2) RayTracer with PyMSES Classic Operator definition #######################
# ####################################################################################
# op = ScalarOperator(lambda dset: dset["rho"])
op = FractionOperator(lambda dset: dset["rho"] ** 2, lambda dset: dset["rho"], ro.info["unit_density"])
# ######## Option A : Create RayTracer = multiprocessing on data #####################
# rt = RayTracer(ro, ["rho"])
# #map = rt.process(op, cam) # multiprocessing cpu on data files loaded
# map = rt.process(op, cam, source=fullOctreeDataSource) # reuse local octree
# ######## Option B : use OctreeRT = multiprocessing on image pixels/rays ############
map, levelmax_map = OctreeRT.process(op, cam, rgb=False)  # reuse OctreeRT from part 1)
# This should give the same result as RayTracer(ro, ["rho"]).process(op, cam)
# mapPylab = apply_log_scale(splatting)
# import pylab as P
# P.imshow(mapPylab)
# P.show()
save_map_HDF5(map, cam, map_name="rt_%s" % outNumber)
save_HDF5_to_img(("./rt_%s.h5" % outNumber), cmap="jet", img_path="./")
os.remove(("./rt_%s.h5" % outNumber))
save_map_HDF5(levelmax_map, cam, map_name="rt_lvlmax_%s" % outNumber)
save_HDF5_to_img(("./rt_lvlmax_%s.h5" % outNumber), cmap="jet", img_path="./", log_sensitive=False)
os.remove(("./rt_lvlmax_%s.h5" % outNumber))
print(("rt total time = %.1f s" % (time() - t0), "mms = ", mms, "max AMR read = ", \
    cam.get_required_resolution()))
