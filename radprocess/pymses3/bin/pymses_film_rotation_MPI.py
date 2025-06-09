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

# basic rotation ray tracing script film
NB_IMG = 180
from time import time
from mpi4py import MPI

comm = MPI.COMM_WORLD
myrank = comm.Get_rank()
total_time = time()
from pymses import RamsesOutput

if myrank == 0: print(("import RamsesOutput=", (time() - total_time)))
import os
from pymses.analysis.visualization import Camera, FractionOperator
from pymses.analysis.visualization.image_plot_utils import *
from pymses.analysis.visualization.raytracing import RayTracerMPI
import numpy as N
from optparse import OptionParser

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
    mms = 256  # 1920
ro = RamsesOutput(fileDir, outNumber)
outNumber = ro.iout
# if myrank == 0:time0 = time()
# if myrank == 0: print "time RamsesOutput=", (time()-time0)
center = [0.5, 0.5, 0.5]
rt = RayTracerMPI(ro, ["rho"], remember_data=True)
# if myrank == 0: print "time RamsesOutput+rt=", (time()-time0)
op = FractionOperator(lambda dset: (dset["rho"] ** 2), lambda dset: (dset["rho"]))
i = 0  # first iteration of the following "for" loop done separately to save the colormap
angle = N.pi * i / NB_IMG + N.pi * .5 / 180  # add .5 to avoid grid alignment artifact problem
axe = [N.sin(angle), 0, N.cos(angle)]
cam = Camera(center=[0.5, 0.5, 0.5], line_of_sight_axis=axe, up_vector="y", region_size=[.5, .28125],
             size_unit=ro.info["unit_length"], distance=0.2, far_cut_depth=0.2, map_max_size=mms)
t0 = time()
map = rt.process(op, cam, use_balanced_cpu_list=True)
t1 = time()
if myrank == 0:
    print(("myrank", myrank, "rt total time = %.1f s" % (t1 - t0), "mms = ", \
        mms, "max AMR read = ", cam.get_required_resolution()))
    save_map_HDF5(map, cam, map_name="img%s" % (i))
    # ran used to fix the colormap during the movie (= use first frame colormap for each frame)
    ran = save_HDF5_to_img("./img%s.h5" % (i), cmap="jet", img_path="./")
    os.remove(("./img%s.h5" % (i)))
for i in range(1, NB_IMG):
    angle = N.pi * i / NB_IMG + N.pi * .5 / 180
    axe = [N.sin(angle), 0, N.cos(angle)]
    cam = Camera(center=[0.5, 0.5, 0.5], line_of_sight_axis=axe, up_vector="y", region_size=[.5, .28125],
                 size_unit=ro.info["unit_length"], distance=0.2, far_cut_depth=0.2, map_max_size=mms)
    t0 = time()
    map = rt.process(op, cam, use_balanced_cpu_list=True)
    t1 = time()
    if myrank == 0:
        print(("myrank", myrank, "rt total time = %.1f s" % (t1 - t0), \
            "mms = ", mms, "max AMR read = ", cam.get_required_resolution()))
        save_map_HDF5(map, cam, map_name="img%s" % (i))
        save_HDF5_to_img("./img%s.h5" % (i), cmap="jet", img_path="./", ran=ran)
        os.remove(("./img%s.h5" % (i)))
if myrank == 0: print(("total film time=", (time() - total_time)))
