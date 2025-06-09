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

# octree multiprocessing ray tracing with local octree movie : 
# galaxy rotation + zoom + rotation + translation
# Local octree reconstruction use : Pay attention to 
# the ngrid_max and misc.NUMBER_OF_PROCESSES_LIMIT parameters
# as you need enough RAM and CPU ressources available !!!

from time import time

myrank = 0
total_time = time()
from pymses import RamsesOutput

if myrank == 0: print "import RamsesOutput=", (time() - total_time)
import os
from pymses.analysis.visualization import Camera, FractionOperator
from pymses.analysis.visualization.image_plot_utils import *
from pymses.analysis.visualization.raytracing import OctreeRayTracer
from pymses.sources.ramses.octree import CameraOctreeDatasource, CameraOctreeDataset
from optparse import OptionParser
import numpy as N

parser = OptionParser()
parser.usage = "%prog mms"
(opts, args) = parser.parse_args()
try:
	mms = int(args[0])
except:
	mms = 192  # 0 # 1920 = full HD movie
from pymses.utils import misc

misc.NUMBER_OF_PROCESSES_LIMIT = 1  # multiprocessing
# This ngrid_max needs to be big enough to fit your data box !
ngrid_max = 2e6  # 10e6 ~ 3.8GB of RAM memory
# If you don't have enough RAM on one node, use an other script with 
# the classic ray tracer or with the MPI ray tracer that can run on many nodes
right_eye = True  # right eye movie
img_start = 0  # Use this to restart a movie computation from this image
img_stop = 10000
nbImgRot = 40
nbImgZoom = 300  # zoom from amr level 13 to 17
nbImgTranslate = 200

# small shift to avoid grid alignment problems
center = [0.5000001, 0.5000001, 0.5000001]
angle = N.pi * .5 / 180  # add .5 to avoid grid alignment artifact problem
los = [0.4, N.sin(angle), N.cos(angle)]
zoom = 2
# following test to see until which level we are going into:
# for i in range(1,nbImg):
#	zoom = zoom*0.99
# print zoom
# cam  = Camera(center=center, line_of_sight_axis=los,
#	region_size=[.25*zoom, .140625*zoom], distance=0.2*zoom, 
#	far_cut_depth=0.2*zoom, map_max_size=mms)
# cam.get_required_resolution()

if myrank == 0: time0 = time()
outnumber = 42
ro = RamsesOutput("/home/labadens", outnumber)
if myrank == 0: print "RamsesOutput time=", (time() - time0)

op = FractionOperator(lambda dset: (dset["rho"] ** 2),
					  lambda dset: (dset["rho"]))
# Colormap initialisation (ran)
cam = Camera(center=center, line_of_sight_axis=los, region_size=[.25 * zoom, .140625 * zoom],
			 size_unit=ro.info["unit_length"], distance=0.2 * zoom, far_cut_depth=0.2 * zoom, map_max_size=mms)

source = ro.amr_source(["rho"])
esize = 0.5 ** (ro.info["levelmin"] + 1)
camOctSource = cam.copy()
fullOctreeDataSource = CameraOctreeDatasource(camOctSource, esize, source, ngrid_max=ngrid_max,
											  include_split_cells=True).dset
OctreeRT = OctreeRayTracer(fullOctreeDataSource)
dataset_already_loaded = True
t0 = time()
map, levelmax_map = OctreeRT.process(op, cam, rgb=False)
t1 = time()
if myrank == 0:
	i, img_number = 0, 0
	print "myrank", myrank, "rt total time = %.1f s" % (t1 - t0), \
		"mms = ", mms, "max AMR read = ", cam.get_required_resolution()
    save_map_HDF5(map, cam, map_name="img%s" % (i))
	# ran used to fix the colormap during the movie
	# (= use first frame colormap for each frame)
	ran = save_HDF5_to_img("./img%s.h5" % (i), cmap="jet", img_path="./")
	os.remove(("./img%s.h5" % (i)))

# ---------------------------------- #
#   Rotate
# ---------------------------------- #
if right_eye:
	istart = 0
else:
	istart = 1  # the first image is already computed for colormap init
for i in range(istart, nbImgRot):
	img_number = i
	compute_img = img_number >= img_start and img_number <= img_stop or \
				  (img_start == 0 and right_eye)
	if compute_img or (i == nbImgRot - 1):
		angle = N.pi * i / nbImgRot * .1 + N.pi * .5 / 180  # add .5 to avoid grid alignment artifact problem
		los = [0.4, N.sin(angle), N.cos(angle)]
		if compute_img:
			cam = Camera(center=center, line_of_sight_axis=los, region_size=[.25 * zoom, .140625 * zoom],
						 size_unit=ro.info["unit_length"], distance=0.2 * zoom, far_cut_depth=0.2 * zoom,
						 map_max_size=mms)
			if right_eye: cam = cam.get_3D_right_eye_cam()
			t0 = time()
			map, levelmax_map = OctreeRT.process(op, cam, rgb=False)
			t1 = time()
			if myrank == 0:
				print "myrank", myrank, "rt total time = %.1f s" % (t1 - t0), \
					"mms = ", mms, "max AMR read = ", cam.get_required_resolution()
				save_map_HDF5(map, cam, map_name="img%s" % (img_number))
				save_HDF5_to_img("./img%s.h5" % (img_number), cmap="jet",
								 img_path="./", ran=ran)
				os.remove(("./img%s.h5" % (img_number)))
# ---------------------------------- #
#   Zoom
# ---------------------------------- #
for i in range(nbImgZoom):
	img_number = i + nbImgRot
	zoom = zoom * 0.99
	compute_img = img_number >= img_start and img_number <= img_stop
	if compute_img or (i == nbImgZoom - 1):
		cam = Camera(center=center, line_of_sight_axis=los, region_size=[.25 * zoom, .140625 * zoom],
					 size_unit=ro.info["unit_length"], distance=0.2 * zoom, far_cut_depth=0.2 * zoom, map_max_size=mms)
		if right_eye: cam = cam.get_3D_right_eye_cam()
		if not camOctSource.contains_camera(cam):
			camOctSource = cam.copy()
			# extend loading by 2 AMR level (anticipate the zoom):
			camOctSource.map_max_size = cam.map_max_size * 4
			fullOctreeDataSource = CameraOctreeDatasource(camOctSource, esize,
														  source, ngrid_max=ngrid_max, include_split_cells=1).dset
			OctreeRT = OctreeRayTracer(fullOctreeDataSource)
			dataset_already_loaded = False
		if compute_img:
			t0 = time()
			map, levelmax_map = OctreeRT.process(op, cam, rgb=False)
			dataset_already_loaded = True
			t1 = time()
			if myrank == 0:
				print "myrank", myrank, "rt total time = %.1f s" % (t1 - t0), "mms = ", \
					mms, "max AMR read = ", cam.get_required_resolution()
				save_map_HDF5(map, cam, map_name="img%s" % (img_number))
				save_HDF5_to_img("./img%s.h5" % (img_number), cmap="jet",
								 img_path="./", ran=ran)
				os.remove(("./img%s.h5" % (img_number)))
# ---------------------------------- #
#   Rotate back
# ---------------------------------- #
for i in range(nbImgRot):
	img_number = i + nbImgRot + nbImgZoom
	compute_img = img_number >= img_start and img_number <= img_stop
	if compute_img or (i == nbImgRot - 1):
		# add .5 to avoid grid alignment artifact problem
		angle = N.pi * (nbImgRot - i - 1) / nbImgRot * .1 + N.pi * .5 / 180
		los = [0.4, N.sin(angle), N.cos(angle)]
		if compute_img:
			cam = Camera(center=center, line_of_sight_axis=los, region_size=[.25 * zoom, .140625 * zoom],
						 size_unit=ro.info["unit_length"], distance=0.2 * zoom, far_cut_depth=0.2 * zoom,
						 map_max_size=mms)
			if right_eye: cam = cam.get_3D_right_eye_cam()
			t0 = time()
			map, levelmax_map = OctreeRT.process(op, cam, rgb=False)
			dataset_already_loaded = True
			t1 = time()
			if myrank == 0:
				print "myrank", myrank, "rt total time = %.1f s" % (t1 - t0), "mms = ", mms, \
					"max AMR read = ", cam.get_required_resolution()
				save_map_HDF5(map, cam, map_name="img%s" % (img_number))
				save_HDF5_to_img("./img%s.h5" % (img_number), cmap="jet", img_path="./", ran=ran)
				os.remove(("./img%s.h5" % (img_number)))
# ---------------------------------- #
#   Translate
# ---------------------------------- #
for i in range(nbImgTranslate):
	img_number = i + 2 * nbImgRot + nbImgZoom
	center[1] += .002 * zoom
	compute_img = img_number >= img_start and img_number <= img_stop
	if compute_img:
		cam = Camera(center=center, line_of_sight_axis=los, region_size=[.25 * zoom, .140625 * zoom],
					 size_unit=ro.info["unit_length"], distance=0.2 * zoom, far_cut_depth=0.2 * zoom, map_max_size=mms)
		if right_eye: cam = cam.get_3D_right_eye_cam()
		if not camOctSource.contains_camera(cam):
			camOctSource = cam.copy()
			# extend the region :
			camOctSource.region_size[1] += .002 * zoom * (nbImgTranslate - i)
			fullOctreeDataSource = CameraOctreeDatasource(camOctSource, esize, source,
														  ngrid_max=ngrid_max, include_split_cells=1).dset
			OctreeRT = OctreeRayTracer(fullOctreeDataSource)
			dataset_already_loaded = False
		t0 = time()
		map, levelmax_map = OctreeRT.process(op, cam, rgb=False)
		dataset_already_loaded = True
		t1 = time()
		if myrank == 0:
			print "myrank", myrank, "rt total time = %.1f s" % (t1 - t0), "mms = ", mms, \
				"max AMR read = ", cam.get_required_resolution()
			save_map_HDF5(map, cam, map_name="img%s" % (img_number))
			save_HDF5_to_img("./img%s.h5" % (img_number), cmap="jet", img_path="./", ran=ran)
			os.remove(("./img%s.h5" % (img_number)))
if myrank == 0: print "total film time=", (time() - total_time)
