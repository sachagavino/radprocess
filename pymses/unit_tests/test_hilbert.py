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
test_hilbert.py -- test module for Hilbert space-filling curve
"""

import pymses.sources.ramses.domain_decomposition as hilbert
import numpy
import os

PWD = os.path.abspath(os.path.dirname(__file__))


def test_hilbert2d():
    points = [[0, 0],
              [1, 0],
              [1, 1],
              [0, 1],
              [0, 2],
              [0, 3],
              [1, 3],
              [1, 2],
              [2, 2],
              [2, 3],
              [3, 3],
              [3, 2],
              [3, 1],
              [2, 1],
              [2, 0],
              [3, 0]]

    keys = numpy.arange(16)

    parray = numpy.array(points)
    hkeys = hilbert.compute_hilbert_key(parray, 4)
    assert (hkeys == keys).all()


def test_hilbert3d():

    # keys of randomly generated cells as computed by the RAMSES hilbert3d
    # routine
    levels = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    points = [[1, 1, 1],
              [2, 2, 2],
              [5, 4, 3],
              [7, 8, 5],
              [9, 5, 2],
              [28, 14, 41],
              [34, 94, 115],
              [96, 128, 113],
              [390, 358, 199],
              [193, 612, 832],
              [1314, 1810, 1444]]

    keys = [5, 46, 286, 1848, 1786, 39597, 617183, 7559591, 76270683,
            326598723, 5623585136]

    for point, level, key in zip(points, levels, keys):
        parray = numpy.array([point])
        hkey = hilbert.compute_hilbert_key(parray, level)
        assert hkey[0] == key
        reverseTest = hilbert.compute_indices(hkey, level, 3)
        assert (reverseTest == parray).all()


def test_positions_to_indices_reverse():
    # keys = numpy.array([104862657821952688, 104862657821953136, 104862657821953472], dtype='uint64')
    # indices = hilbert.compute_indices(keys, 19, 3)
    indices = numpy.array([[297695, 307258, 293165],
                           [297691, 307258, 293157],
                           [297688, 307260, 293152]])
    order_list = numpy.ones(3, 'i') * 19
    positions = hilbert.indices_to_positions(indices, order_list)
    indices_recovered = hilbert.positions_to_indices(positions, 19)
    assert (indices_recovered == indices).all()


# point_keys = hilbert.compute_hilbert_key(indices, 19)
# assert (point_keys == keys).all()

# def data_test_AQUARIUS_DATA_map_points_error_count():
#	import pymses
#	ro = pymses.RamsesOutput("/data/Aquarius/output",193)
#	source = ro.amr_source(["rho"])
#	icpuTest = 50
#	dset = source.get_domain_dset(icpuTest)
#	pts_act = dset.get_cell_centers()[dset.get_active_mask()]
#	pts_act = (pts_act.reshape((len(pts_act)*8,3)))
#	dom = source.dom_decomp
#	map_points = dom.map_points(pts_act)
#	error_count = sum(map_points!=icpuTest)
#	print "errors counted =", error_count, " out of", len(pts_act), "points"
#	assert (error_count == 0)


def oldtest_map_box_512():

    # Obtained by running amr2cube
    keys_min, keys_max = numpy.loadtxt(
        os.path.join(PWD, "hilbert_keys_test_512.dat"), unpack=True)
    levelmax = 17
    bbox_min = [0.2, 0.3, 0.6]
    bbox_max = [0.23, 0.45, 0.67]
    cpu_list = [92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105,
                106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118,
                119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131,
                132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144,
                145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157,
                158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170,
                171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183,
                184, 185]

    hdcp = hilbert.HilbertDomainDecomp(3, keys_min, keys_max, (0, levelmax))
    pymses_cpu_list = hdcp.map_box((bbox_min, bbox_max))
    print(pymses_cpu_list)
    print(cpu_list)
    assert numpy.all(pymses_cpu_list == numpy.array(cpu_list, 'i'))


def test_map_box_32():

    # Obtained by running amr2cube
    keys_min, keys_max = numpy.loadtxt(
        os.path.join(PWD, "hilbert_keys_test_32.dat"), unpack=True)
    levelmax = 13
    bbox_min = 3 * [0.45]
    bbox_max = 3 * [0.60]
    cpu_list = [1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 21,
                22, 23, 24, 25, 26, 27, 29, 30, 31, 32]

    hdcp = hilbert.HilbertDomainDecomp(3, keys_min, keys_max, (0, levelmax))
    pymses_cpu_list = hdcp.map_box((bbox_min, bbox_max))
    assert numpy.all(pymses_cpu_list == numpy.array(cpu_list, 'i'))


def test_minimal_domain():
    # Domain decomposition 5 process (Red, Violet, Green, Orange, Blue) on a 16ix16 grid
    # -----------------------------------------------------------------------------------

    ###########################################################################
    #    0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16    #
    # 16 .---------------------------------------------------------------. 16 #
    #    | G | G | G | G | G | G | G | G | G | G | O | O | O | O | O | O |    #
    # 15 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---| 15 #
    #    | G | G | G | G | G | G | G | G | G | O | O | O | O | O | O | O |    #
    # 14 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---| 14 #
    #    | G | G | G | G | G | G | G | G | G | G | O | O | O | O | O | O |    #
    # 13 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---| 13 #
    #    | V | G | G | G | G | G | G | G | G | G | O | O | O | O | O | O |    #
    # 12 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---| 12 #
    #    | V | V | V | V | G | G | G | G | G | G | G | G | O | O | O | O |    #
    # 11 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---| 11 #
    #    | V | V | V | V | G | G | G | G | G | G | G | G | O | O | O | O |    #
    # 10 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---| 10 #
    #    | V | V | V | V | G | G | G | G | G | G | G | G | O | O | O | O |    #
    #  9 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|  9 #
    #    | V | V | V | V | G | G | G | G | G | G | G | G | O | O | O | O |    #
    #  8 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|  8 #
    #    | V | V | V | V | V | V | V | V | O | O | O | O | O | O | O | O |    #
    #  7 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|  7 #
    #    | V | V | V | V | V | V | V | V | O | O | O | O | O | O | O | O |    #
    #  6 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|  6 #
    #    | V | V | V | V | R | R | V | V | O | O | B | B | O | O | O | O |    #
    #  5 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|  5 #
    #    | V | V | V | V | R | R | V | V | O | O | B | B | O | O | O | O |    #
    #  4 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|  4 #
    #    | R | R | R | R | R | R | R | R | B | B | B | B | B | B | B | B |    #
    #  3 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|  3 #
    #    | R | R | R | R | R | R | R | R | B | B | B | B | B | B | B | B |    #
    #  2 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|  2 #
    #    | R | R | R | R | R | R | R | R | B | B | B | B | B | B | B | B |    #
    #  1 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|  1 #
    #    | R | R | R | R | R | R | R | R | B | B | B | B | B | B | B | B |    #
    #  0 '---------------------------------------------------------------'  0 #
    #    0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16    #
    ###########################################################################
    ndim = 2
    cell_keymin = numpy.array([0.0, 36.0, 81.0, 151.0, 220.0])
    cell_keymax = numpy.array([36.0, 81.0, 151.0, 220.0, 256.0])
    nlevelmax = 4

    corner_keymin = cell_keymin * 2 ** ndim
    corner_keymax = cell_keymax * 2 ** ndim
    dd = hilbert.HilbertDomainDecomp(ndim, corner_keymin, corner_keymax, (1, nlevelmax))

    # Domain 0 #
    idomain = 0
    grids, orders, no = dd.minimal_domain(idomain, read_lmax=3)
    assert (grids * 2 ** nlevelmax == numpy.array([[2., 2.],
                                                   [6., 2.],
                                                   [5., 5.]])).all()
    assert (2. ** (4 - orders) == numpy.array([4., 4., 2.])).all()
    assert no == 0

    # Domain 1 #
    idomain = 1
    grids, orders, no = dd.minimal_domain(idomain, read_lmax=3)
    assert (grids * 2 ** nlevelmax == numpy.array([[2., 6.],
                                                   [2., 10.],
                                                   [7., 5.],
                                                   [5., 7.],
                                                   [7., 7.],
                                                   [1., 13.]])).all()
    assert (2. ** (4 - orders) == numpy.array([4., 4., 2., 2., 2., 2.])).all()
    assert no == 1

    # Domain 2 #
    idomain = 2
    grids, orders, no = dd.minimal_domain(idomain, read_lmax=3)
    assert (grids * 2 ** nlevelmax == numpy.array([[6., 10.],
                                                   [6., 14.],
                                                   [10., 10.],
                                                   [3., 13.],
                                                   [1., 15.],
                                                   [3., 15.],
                                                   [9., 13.],
                                                   [1., 13.],
                                                   [9., 15.]])).all()
    assert (2. ** (4 - orders) == numpy.array([4., 4., 4., 2., 2., 2., 2., 2., 2.])).all()
    assert no == 2

    # Domain 3 #
    idomain = 3
    grids, orders, no = dd.minimal_domain(idomain, read_lmax=3)
    assert (grids * 2 ** nlevelmax == numpy.array([[14., 6.],
                                                   [14., 10.],
                                                   [14., 14.],
                                                   [9., 5.],
                                                   [9., 7.],
                                                   [11., 7.],
                                                   [11., 13.],
                                                   [11., 15.],
                                                   [9., 15.]])).all()
    assert (2. ** (4 - orders) == numpy.array([4., 4., 4., 2., 2., 2., 2., 2., 2.])).all()
    assert no == 1

    # Domain 4 #
    idomain = 4
    grids, orders, no = dd.minimal_domain(idomain, read_lmax=3)
    assert (grids * 2 ** nlevelmax == numpy.array([[10., 2.],
                                                   [14., 2.],
                                                   [11., 5.]])).all()
    assert (2. ** (4 - orders) == numpy.array([4., 4., 2.])).all()
    assert no == 0


def test_map_box():
    # Domain decomposition 5 process (Red, Violet, Green, Orange, Blue) on a 16ix16 grid
    # -----------------------------------------------------------------------------------

    ###########################################################################
    #    0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16    #
    # 16 .---------------------------------------------------------------. 16 #
    #    | G | G | G | G | G | G | G | G | G | G | O | O | O | O | O | O |    #
    # 15 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---| 15 #
    #    | G | G | G | G | G | G | G | G | G | O | O | O | O | O | O | O |    #
    # 14 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---| 14 #
    #    | G | G | G | G | G | G | G | G | G | G | O | O | O | O | O | O |    #
    # 13 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---| 13 #
    #    | V | G | G | G | G | G | G | G | G | G | O | O | O | O | O | O |    #
    # 12 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---| 12 #
    #    | V | V | V | V | G | G | G | G | G | G | G | G | O | O | O | O |    #
    # 11 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---| 11 #
    #    | V | V | V | V | G | G | G | G | G | G | G | G | O | O | O | O |    #
    # 10 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---| 10 #
    #    | V | V | V | V | G | G | G | G | G | G | G | G | O | O | O | O |    #
    #  9 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|  9 #
    #    | V | V | V | V | G | G | G | G | G | G | G | G | O | O | O | O |    #
    #  8 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|  8 #
    #    | V | V | V | V | V | V | V | V | O | O | O | O | O | O | O | O |    #
    #  7 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|  7 #
    #    | V | V | V | V | V | V | V | V | O | O | O | O | O | O | O | O |    #
    #  6 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|  6 #
    #    | V | V | V | V | R | R | V | V | O | O | B | B | O | O | O | O |    #
    #  5 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|  5 #
    #    | V | V | V | V | R | R | V | V | O | O | B | B | O | O | O | O |    #
    #  4 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|  4 #
    #    | R | R | R | R | R | R | R | R | B | B | B | B | B | B | B | B |    #
    #  3 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|  3 #
    #    | R | R | R | R | R | R | R | R | B | B | B | B | B | B | B | B |    #
    #  2 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|  2 #
    #    | R | R | R | R | R | R | R | R | B | B | B | B | B | B | B | B |    #
    #  1 |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|  1 #
    #    | R | R | R | R | R | R | R | R | B | B | B | B | B | B | B | B |    #
    #  0 '---------------------------------------------------------------'  0 #
    #    0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16    #
    ###########################################################################
    ndim = 2
    cell_keymin = numpy.array([0.0, 36.0, 81.0, 151.0, 220.0])
    cell_keymax = numpy.array([36.0, 81.0, 151.0, 220.0, 256.0])
    nlevelmax = 4

    corner_keymin = cell_keymin * 2 ** ndim
    corner_keymax = cell_keymax * 2 ** ndim
    dd = hilbert.HilbertDomainDecomp(ndim, corner_keymin, corner_keymax, (1, nlevelmax))

    bbox_min = [0., 0.]
    bbox_max = [0.25, 0.25]
    cpu_list = [1]

    pymses_cpu_list = dd.map_box((bbox_min, bbox_max))
    assert numpy.all(pymses_cpu_list == numpy.array(cpu_list, 'i'))

    bbox_min = [0., 0.]
    bbox_max = [0.5, 0.5]
    cpu_list = [1, 2]

    pymses_cpu_list = dd.map_box((bbox_min, bbox_max))
    assert numpy.all(pymses_cpu_list == numpy.array(cpu_list, 'i'))

    bbox_min = [0.375, 0.375]
    bbox_max = [0.625, 0.625]
    cpu_list = [2, 3, 4]

    pymses_cpu_list = dd.map_box((bbox_min, bbox_max))
    assert numpy.all(pymses_cpu_list == numpy.array(cpu_list, 'i'))
