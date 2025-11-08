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


import numpy as N
import os
from pymses.analysis import Camera
from pymses.utils import constants as C
from nose import tools


class TestCameraHDF5:
    _h5f = None

    @classmethod
    def setup_class(cls):
        import h5py
        cls._pwd = os.path.dirname(__file__)
        cls._h5_dir = os.path.join(cls._pwd, "pymses4_old")
        cls._h5fname_legacy = os.path.join(cls._h5_dir, "Aquarius_density_RT_3disks.h5")
        cls._h5_dir_v1 = os.path.join(cls._pwd, "camera_v1")
        cls._h5F_v1 = h5py.File(os.path.join(cls._h5_dir_v1, "camera.h5"), 'r')

    def test_old_3disks_camera(self):
        print("Test old format camera HDF5 input")
        h5cam = Camera.load_legacy_HDF5_camera(self._h5fname_legacy)
        center = [0.56820763392224549, 0.5868017758212859, 0.56033614395861875]
        los = [0.12682409771262307, -0.98100954931698847, -0.14675119211869589]
        vaxis = [0.30198897499497745, -0.10273635375054604, 0.94775941071535974]
        bounds = N.array([[-0.0022674599999999981, -0.0022674599999999981, -0.0022674599999999981],
                          [0.0022674599999999981, 0.0022674599999999981, 0.0022674599999999981]])
        fcd = -bounds[0, 2]
        dist = bounds[1, 2]
        rsize = N.diff(bounds, axis=0)[0, :-1]

        cam = Camera(center=center, line_of_sight_axis=los, up_vector=vaxis, region_size=rsize, size_unit=C.pc,
                     distance=dist, far_cut_depth=fcd, map_max_size=1030, log_sensitive=True)
        tools.assert_equals(h5cam, cam)

    def test_camera_v1_hdf5_io(self):
        print("Test Camera object HDF5 I/O")
        center = [0.41254, 0.60254, 0.74542]
        los = [0.3019889750, -0.10273635375, 0.94775941]
        vaxis = [0.1268240977, -0.981009549317, -0.146751192]

        # C.Rsun has latex format string
        cam = Camera(center=center, line_of_sight_axis=los, up_vector=vaxis, region_size=[0.3, 0.45], size_unit=C.Rsun,
                     distance=0.3, far_cut_depth=0.25, map_max_size=600, log_sensitive=True)

        # cam.save_HDF5(self._h5F_v1)
        cam2 = Camera.from_HDF5(self._h5F_v1)
        tools.assert_equal(cam, cam2)

    @classmethod
    def teardown_class(cls):
        if cls._h5F_v1 is not None:
            cls._h5F_v1.close()
