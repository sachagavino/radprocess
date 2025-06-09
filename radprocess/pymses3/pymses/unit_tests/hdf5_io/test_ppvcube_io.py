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
from nose import tools
import tempfile

from pymses import RamsesOutput
from pymses.analysis.ppv import PPVCube, PPVProcessor
from pymses.analysis import Camera


class TestPPVCubeHDF5:
    @classmethod
    def setup_class(cls):
        cls._pwd = os.path.dirname(__file__)
        cls._h5_dir_v1 = os.path.join(cls._pwd, "ppvcube_v1")
        cls._h5fname_v1 = os.path.join(cls._h5_dir_v1, "Aquarius_ppv_v1.h5")
        fd, cls._h5f_v1_temp = tempfile.mkstemp(prefix="pymses_ppvcube", suffix=".h5")
        os.close(fd)
        cls._ro = RamsesOutput(os.path.join(os.path.expanduser("~"), "data/Aquarius/output"), 193)

    def test_ppvcube_v1_hdf5_io(self):
        print("Test PPVCube object HDF5 I/O")
        processor = PPVProcessor(self._ro, "rho", "vel")
        cam = Camera(center=[0.566561, 0.583253, 0.555024], map_max_size=64, region_size=[0.023, 0.023],
                     distance=0.0115, far_cut_depth=0.0115)
        cube = processor.process(cam, vrange=[-450.0, -50.0], nv=20)
        # Save computed PPVCube object into temporary HDF5 file
        cube.save_HDF5(self._h5f_v1_temp, float32=True)
        del cube

        # Reload PPVCube object from temporary HDF5 file
        cube_temp = PPVCube.from_HDF5(self._h5f_v1_temp)

        # Compare it to the reference PPVCube object
        cube_true = PPVCube.from_HDF5(self._h5fname_v1)
        tools.assert_equal(cube_true, cube_temp)

    @classmethod
    def teardown_class(cls):
        if os.path.isfile(cls._h5f_v1_temp):
            os.remove(cls._h5f_v1_temp)
