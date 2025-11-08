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
import matplotlib.pyplot as P
import os

from pymses import RamsesOutput
from pymses.analysis.ppv.processor import PPVProcessor
from pymses.analysis import Camera


class TestPPVRayCast(object):
    _datadir = None
    _iout = None
    _ro = None

    @classmethod
    def setup_class(cls):
        cls._datadir = os.path.join(os.path.expanduser("~"), "data/Aquarius/output")
        cls._iout = 193
        cls._ro = RamsesOutput(cls._datadir, cls._iout)

    def test_ppv_cube(self):
        processor = PPVProcessor(self._ro, "rho", "vel")
        cam = Camera(map_max_size=1024, region_size=[0.25, 0.25])
        cube = processor.process(cam, vrange=[-100.0, 100.0], nv=50)
        v_hist = N.array([2.66119128e+17, 7.44225948e+17, 1.18754462e+17, 9.08913705e+17])
        tools.assert_true(N.allclose(cube.data[512, 512, 4:8], v_hist, rtol=1.0e-6))
        # cube.save_fits("my_ppv")

    @classmethod
    def teardown_class(cls):
        pass
