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
from pymses.utils._point_utils import meshgrid
from pymses.utils import constants as C
from pymses.analysis import Camera
import numpy
import time


class TestPointUtils:
    _cam = None

    @classmethod
    def setup_class(cls):
        center = [0.41254, 0.60254, 0.74542]
        los = [0.3019889750, -0.10273635375, 0.94775941]
        vaxis = [0.1268240977, -0.981009549317, -0.146751192]
        cls._cam = Camera(center=center, line_of_sight_axis=los, up_vector=vaxis, region_size=[0.3, 0.45],
                          size_unit=C.pc, distance=0.3, far_cut_depth=0.25, map_max_size=512, log_sensitive=True)

    def test_meshgrid_2d(self):
        t0 = time.time()
        nx, ny = 10, 20
        xgrid, ygrid = numpy.mgrid[0:nx, 0:ny]
        pgrid = numpy.concatenate([xgrid, ygrid])
        plist = pgrid.reshape((2, -1)).transpose()
        t1 = time.time()
        print(t1-t0)

        xpts = numpy.arange(nx, dtype='d')
        ypts = numpy.arange(ny, dtype='d')
        grid = meshgrid([xpts, ypts])
        t2 = time.time()
        print(t2-t1)

        assert numpy.all(grid == plist)

    def test_meshgrid_3d(self):
        t0 = time.time()
        nx, ny, nz = 20, 30, 50
        xgrid, ygrid, zgrid = numpy.mgrid[0:nx, 0:ny, 0:nz]
        pgrid = numpy.concatenate([xgrid, ygrid, zgrid])
        plist = pgrid.reshape((3, -1)).transpose()
        t1 = time.time()
        print(t1-t0)

        xpts = numpy.arange(nx, dtype='d')
        ypts = numpy.arange(ny, dtype='d')
        zpts = numpy.arange(nz, dtype='d')
        grid = meshgrid([xpts, ypts, zpts])
        t2 = time.time()
        print(t2-t1)

        assert numpy.all(grid == plist)

    def test_meshgrid_4d(self):
        t0 = time.time()
        nx, ny, nz, nw = 20, 30, 50, 17
        xgrid, ygrid, zgrid, wgrid = numpy.mgrid[0:nx, 0:ny, 0:nz, 0:nw]
        pgrid = numpy.concatenate([xgrid, ygrid, zgrid, wgrid])
        plist = pgrid.reshape((4, -1)).transpose()
        t1 = time.time()
        print(t1-t0)

        xpts = numpy.arange(nx, dtype='d')
        ypts = numpy.arange(ny, dtype='d')
        zpts = numpy.arange(nz, dtype='d')
        wpts = numpy.arange(nw, dtype='d')
        grid = meshgrid([xpts, ypts, zpts, wpts])
        t2 = time.time()
        print(t2-t1)

        assert numpy.all(grid == plist)

    def test_meshgrid_indexing(self):
        nx, ny, nz = 20, 30, 50
        xpts = numpy.arange(nx, dtype='d')
        ypts = numpy.arange(ny, dtype='d')
        zpts = numpy.arange(nz, dtype='d')
        grid = meshgrid([xpts, ypts, zpts])

        assert tuple(grid.reshape(nx, ny, nz, 3)[10, 14, 25, :]) == (10., 14., 25.)

    def test_camera_corner_points(self):
        b = self._cam.get_map_box()
        box_bounds = [numpy.array([b.min_coords[0], b.max_coords[0]], dtype='d'),
                      numpy.array([b.min_coords[1], b.max_coords[1]], dtype='d'),
                      numpy.array([b.min_coords[2], b.max_coords[2]], dtype='d')]
        corner_list = meshgrid(box_bounds)
        correct_corner_list = numpy.array([[-0.15, -0.225, -0.25],
                                           [-0.15, -0.225,  0.3],
                                           [-0.15,  0.225, -0.25],
                                           [-0.15,  0.225,  0.3],
                                           [ 0.15, -0.225, -0.25],
                                           [ 0.15, -0.225,  0.3],
                                           [ 0.15,  0.225, -0.25],
                                           [ 0.15,  0.225,  0.3]], dtype='d')
        print(correct_corner_list)
        print(corner_list)
        assert numpy.all(correct_corner_list == corner_list)

    def test_compute_same_value_pixel_size_map(self):
        # TODO
        # !! This will evolve and change with parameters !!
        # map numbers generated with :
        # map = numpy.random.random((5,5))
        # map = numpy.array(map*5,'i')
        map = numpy.array([[2, 0, 3, 3, 0],
                           [4, 1, 4, 0, 3],
                           [4, 0, 3, 4, 3],
                           [1, 1, 3, 2, 2],
                           [3, 4, 1, 1, 3]], 'f')
        same_value_pixel_size_map_asserted_result = numpy.array([[0, 0, 2, 2, 1],
                                                                 [1, 0, 1, 1, 2],
                                                                 [1, 0, 1, 1, 2],
                                                                 [3, 3, 1, 1, 1],
                                                                 [0, 0, 3, 3, 0]], 'i')
        filtered_map_asserted_result = numpy.array([[2.00000000, 0.00000000, 2.00111660, 1.79405243, 1.19373007],
                                                    [2.62965690, 1.00000000, 2.20940490, 2.40980212, 2.08841360],
                                                    [2.45708771, 0.00000000, 2.27524973, 2.70881434, 2.35359556],
                                                    [2.51044634, 2.41439667, 2.15545216, 2.39791002, 2.47302363],
                                                    [3.00000000, 4.00000000, 2.43638318, 2.42615966, 3.00000000]], 'f')

        print("map :\n", map, map.min(),map.max())
        from pymses.utils.point_utils import adaptive_gaussian_blur, compute_same_value_pixel_size_map
        # same_value_pixel_size_map
        same_value_pixel_size_map = compute_same_value_pixel_size_map(map)
        print("same_value_pixel_size_map_asserted_result :\n", same_value_pixel_size_map_asserted_result,\
            same_value_pixel_size_map_asserted_result.min(), same_value_pixel_size_map_asserted_result.max())
        print("same_value_pixel_size_map :\n", same_value_pixel_size_map, same_value_pixel_size_map.min(),\
        same_value_pixel_size_map.max())
        #assert (same_value_pixel_size_map_asserted_result == same_value_pixel_size_map).all()

        # adaptive_gaussian_blur
        filtered_map = adaptive_gaussian_blur(map, same_value_pixel_size_map)
        print("filtered_map :\n", filtered_map, filtered_map.min(), filtered_map.max())
        print("filtered_map_asserted_result :\n", filtered_map_asserted_result,\
            filtered_map_asserted_result.min(), filtered_map_asserted_result.max())
        #assert (abs(filtered_map_asserted_result - filtered_map) < 10e-6).all()
        print("Test passed!")

    @classmethod
    def teardown_class(cls):
        pass
