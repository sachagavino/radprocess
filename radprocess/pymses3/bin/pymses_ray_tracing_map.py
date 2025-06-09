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

# ray tracing map script
import pymses, os
from pymses.analysis.visualization import *
from pymses.analysis.visualization.image_plot_utils import *
from time import time
from pymses.analysis.visualization.raytracing import RayTracer, RayTracerMPI
from optparse import OptionParser

t0 = time()
parser = OptionParser()
parser.usage = "%prog ramses_directory ramses_output_number map_max_resolution=512"
(opts, args) = parser.parse_args()
try:
    fileDir = args[0]
    outNumber = int(args[1])
except:
    # "None" leads to an automatic look for
    # RAMSES output in the current directory
    fileDir = None
    outNumber = None
try:
    mms = int(args[2])
except:
    mms = 512
ro = pymses.RamsesOutput(fileDir, outNumber)
outNumber = ro.iout
cam = Camera(center=[0.5, 0.5, 0.5], line_of_sight_axis="z", up_vector="y", region_size=[5.0E-1, 5.0E-1],
             size_unit=ro.info["unit_length"], distance=2.5E-1, far_cut_depth=2.5E-1, map_max_size=mms)
# import tables
# file= tables.openFile("camera2Zoom314.h5", "r")
# cam = Camera.from_HDF5(file)
# file.close()
op = FractionOperator(lambda dset: (dset["rho"] ** 2), lambda dset: (dset["rho"]))
# op = ScalarOperator(lambda dset: dset["rho"])
# cam.log_sensitive=False
# op = MaxLevelOperator()
try:  # Use MPI if possible:
    useMPI = False
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    MPI_process_number = comm.Get_size()
    if MPI_process_number == 1:
        print("Try \"mpirun -n 4 pymses_ray_tracing_map.py\" to use 4 MPI raytracing process")
        raise Exception()
    else:
        useMPI = True
except:
    myrank = 0
    rt = RayTracer(ro, ["rho"])
    map = rt.process(op, cam)
    print("Currently using multiprocessing rayTracer")
if useMPI:
    if myrank == 0: print(("Import time =", time() - t0, "MPI_process_number = ", MPI_process_number))
    rt = RayTracerMPI(ro, ["rho"])
    map = rt.process(op, cam)  # , use_balanced_cpu_list=True)
if myrank == 0:
    t1 = time()
    print(("rt total time = %.1f s" % (t1 - t0), "mms = ", mms, "max AMR read = ", cam.get_required_resolution()))
    # mapPylab = apply_log_scale(map)
    # import pylab as P
    # P.imshow(mapPylab)
    # P.show()
    save_map_HDF5(map, cam, map_name="rt_%s" % outNumber)
    save_HDF5_to_img(("./rt_%s.h5" % outNumber), cmap="jet", img_path="./", ramses_output=ro)
# os.remove(("./rt_%s.h5"%outNumber))
