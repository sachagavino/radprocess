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
"""
test_octree_dset.py -- test module for Hilbert space-filling curve
"""

import pymses
from pymses.analysis.splatting import SplatterProcessor
from pymses.analysis.raytracing import RayTracer
# from pymses.analysis.transfer_functions import ColorLinesTransferFunction
from pymses.analysis import Camera, ScalarOperator, slicing, raytracing
from pymses.utils import constants as C
import numpy
from nose import tools
import os


class TestOctreeDataset:
    _h5f = None
    _h5f_v1 = None

    @classmethod
    def setup_class(cls):
        import h5py
        PWD = os.path.abspath(os.path.dirname(__file__))
        cls._octree_amr_dset_h5f = h5py.File(os.path.join(PWD, "test_amr_dset.h5"), 'r')
        cls.asserted_SliceMap_result = numpy.array([[-5.98442495, -5.98442495, -6.00830233],
                                                    [-5.98442495, -5.98442495, -6.00830233],
                                                    [-6.03202943, -6.03202943, -6.07758903]])
        cls.asserted_interpolated_SliceMap_result = numpy.array([[-5.99014809, -5.96351829, -5.96968206],
                                                                 [-6.01119872, -5.99727729, -5.99045202],
                                                                 [-6.02887729, -6.01696159, -6.01550138]])
        cls.asserted_rt_map = numpy.array([[-4.61547079, -4.50387064, -4.50387064],
                                           [-5.17418872, -4.88956574, -4.88956574],
                                           [-5.17418872, -4.88956574, -4.88956574]])

        cls.asserted_rotated_rt_map = numpy.array([[4.43912944e-06, 9.40160105e-06, 6.78405933e-06],
                                                   [5.16985554e-06, 3.88652536e-05, 4.32777455e-05],
                                                   [2.06883743e-06, 7.71202140e-06, 8.58699623e-06]])

        cls.asserted_lvlmax_map = numpy.array([[3, 3, 3],
                                               [3, 3, 3],
                                               [3, 3, 3]])

        cls.asserted_splatted_map = numpy.array([[-0.22320506, -0.0246566, -0.11069895],
                                                 [-0.2586495, -0.0441908, -0.11421175],
                                                 [-0.67037297, -0.41799685, -0.45085401]])

        cls.ramses_output_info = {"levelmin": 3, "levelmax": 3, "dom_decomp": None}

    def test_loading_OctreeDataset(self):
        print("Test AMR octree dataset loading")
        amr_source = pymses.sources.ramses.octree.OctreeDataset.from_HDF5(self._octree_amr_dset_h5f)
        tools.assert_equal(amr_source.get_read_levelmax(), 3)
        tools.assert_equal(amr_source.data_list, [0])
        tools.assert_equal(amr_source.fields["rho"].shape[0], 73)
        tools.assert_true(numpy.allclose(amr_source["rho"][5, 5], 1.33670771473e-07, atol=1.0e-6))

    def test_SliceMap(self):
        print("Test slice map over an AMR octree dataset")
        amr_source = pymses.sources.ramses.octree.OctreeDataset.from_HDF5(self._octree_amr_dset_h5f)
        rho_op = ScalarOperator(lambda dset: dset["rho"], C.H_cc)
        cam = Camera(center=[0.5, 0.5, 0.5], line_of_sight_axis='z', up_vector='y', region_size=[.3, .3], distance=0.,
                     far_cut_depth=0., map_max_size=3, log_sensitive=True)
        datamap = slicing.SliceMap(amr_source, cam, rho_op, z=0.3)
        tools.assert_true(numpy.allclose(numpy.log10(datamap.map), self.asserted_SliceMap_result, atol=1e-6))

    def test_SliceMap_interpolation(self):
        print("Test interpolated slice map over an AMR octree dataset")
        amr_source = pymses.sources.ramses.octree.OctreeDataset.from_HDF5(self._octree_amr_dset_h5f)
        rho_op = ScalarOperator(lambda dset: dset["rho"], C.H_cc)
        cam = Camera(center=[0.5, 0.5, 0.5], line_of_sight_axis='z', region_size=[.3, .3], distance=0.,
                     far_cut_depth=0., up_vector='y', map_max_size=3, log_sensitive=True)
        datamap = slicing.SliceMap(amr_source, cam, rho_op, z=0.3, interpolation=True)
        tools.assert_true(numpy.allclose(numpy.log10(datamap.map), self.asserted_interpolated_SliceMap_result,
                                         atol=1e-6))

    def test_RayTracer(self):
        print("Test ray-traced map over an AMR octree dataset")
        octree_dset = pymses.sources.ramses.octree.OctreeDataset.from_HDF5(self._octree_amr_dset_h5f)
        rho_op = ScalarOperator(lambda dset: dset["rho"], C.H_cc)
        cam = Camera(center=[0.5, 0.5, 0.5], line_of_sight_axis='z', up_vector='y', region_size=[.3, .3], distance=0.3,
                     far_cut_depth=0.3, map_max_size=3, log_sensitive=False)

        rt = raytracing.RayTracer(octree_dset, self.ramses_output_info, rho_op)
        # rt.disable_multiprocessing()
        dmap = rt.process(cam, surf_qty=False)
        # print dmap.map, self.asserted_rt_map
        tools.assert_true(numpy.allclose(numpy.log10(dmap.map), self.asserted_rt_map, atol=1e-6))

    def test_RayTracer_rotated(self):
        print("Test rotated ray-traced map over an AMR octree dataset")
        octree_dset = pymses.sources.ramses.octree.OctreeDataset.from_HDF5(self._octree_amr_dset_h5f)
        rho_op = ScalarOperator(lambda dset: dset["rho"], C.H_cc)
        # camera center small shift is still with OctreeRayTracer to avoid grid limit pb
        cam = Camera(center=[0.501, 0.501, 0.501], line_of_sight_axis=[0.401, 0.601, -0.701], up_vector='y',
                     region_size=[.3, .3], distance=0.3, far_cut_depth=0.3, map_max_size=3, log_sensitive=False)

        rt = raytracing.RayTracer(octree_dset, self.ramses_output_info, rho_op)
        dmap = rt.process(cam, surf_qty=False)
        tools.assert_true(numpy.allclose(dmap.map, self.asserted_rotated_rt_map, atol=1e-3))

    def test_Splatting(self):
        print("Test splatted map over an AMR octree dataset")
        octree_dset = pymses.sources.ramses.octree.OctreeDataset.from_HDF5(self._octree_amr_dset_h5f)
        rho_op = ScalarOperator(lambda dset: dset["rho"], C.H_cc)
        cam = Camera(center=[0.5, 0.5, 0.5], line_of_sight_axis='z', up_vector='y', region_size=[.3, .3], distance=0.3,
                     far_cut_depth=0.3, map_max_size=3, log_sensitive=True)
        mp = SplatterProcessor(octree_dset, self.ramses_output_info, rho_op)
        # mp.disable_multiprocessing()
        map = mp.process(cam, pre_flatten=False, surf_qty=False)

        # map = ImgPlot.apply_log_scale(map)
        # print numpy.log10(map.map)
        # print self.asserted_map
        tools.assert_true(numpy.allclose(numpy.log10(map.map), self.asserted_splatted_map, atol=1e-6))

    @classmethod
    def teardown_class(cls):
        if cls._octree_amr_dset_h5f is not None:
            cls._octree_amr_dset_h5f.close()


__all__ = ['TestOctreeDataset']
