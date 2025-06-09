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

from nose import tools
import numpy as N
import os
import tempfile
from PIL import Image as I

from pymses import RamsesOutput
from pymses.analysis import MassWeightedDensityOperator, raytracing, SphericalCamera, DataMap
from pymses.utils import constants as C


class TestSphericalRayCast(object):
    _datadir = None
    _iout = None
    _ro = None

    @classmethod
    def setup_class(cls):
        cls._datadir = os.path.join(os.path.expanduser("~"), "data/MC_RT2")
        cls._iout = 240
        cls._ro = RamsesOutput(cls._datadir, cls._iout)
        cls._pwd = os.path.dirname(__file__)
        cls._h5_dir_v1 = os.path.join(cls._pwd, "hdf5_io", "datamap_v1")
        cls._h5f_v1 = os.path.join(cls._h5_dir_v1, "MC_RT2_00240_spherical_1024x512.h5")
        cls._v1_img = os.path.join(cls._h5_dir_v1, "MC_RT2_00240_spherical_1024x512_thermal.png")
        fd, cls._h5f_v1_temp = tempfile.mkstemp(prefix="pymses_datamap", suffix=".h5")
        os.close(fd)

    def test_spherical_map(self):
        # Camera
        c = [0.504, 0.506, 0.747]
        los = [-0.5, 0.0, 0.866]
        up = [0.866, 0.0, 0.5]
        nazimuth = 1024
        cam = SphericalCamera(center=c, line_of_sight_axis=los, up_vector=up, region_size=[2.0, 2.0],
                              distance=0.8, far_cut_depth=0.8,
                              log_sensitive=True, map_max_size=nazimuth, size_unit=30.349659078895829 * C.pc)
        lcam, rcam = cam.get_3D_cameras(perspective_angle_deg=2.0)

        amr = self._ro.amr_source(["rho"])

        # Map operator : gas mass-weighted density
        scal_func = MassWeightedDensityOperator("rho", self._ro.info["unit_density"])

        # Map processing
        rt = raytracing.RayTracer(amr, self._ro.info, scal_func)
        # rt.disable_multiprocessing()
        ldmap = rt.process(lcam)

        ldmap.save_HDF5(self._h5f_v1_temp)
        loaded_ldmap = DataMap.from_HDF5(self._h5f_v1_temp)
        ref_ldmap = DataMap.from_HDF5(self._h5f_v1)

        # DataMap comparison
        tools.assert_equal(loaded_ldmap, ldmap)
        tools.assert_equal(loaded_ldmap, ref_ldmap)

        ref_img = I.open(self._v1_img, 'r')
        save_img = loaded_ldmap.save_PNG(cmap="Dark_Thermal")

        a_refimg = N.array(ref_img.getdata())
        a_simg = N.array(save_img.getdata())
        tools.assert_true(N.allclose(a_refimg, a_simg, rtol=1.0e-4))

    @classmethod
    def teardown_class(cls):
        if os.path.isfile(cls._h5f_v1_temp):
            os.remove(cls._h5f_v1_temp)
