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
import numpy
from pymses.utils.misc import *


def test_concatenate_reorder():

    alist = [numpy.random.uniform(size=[3, 2, 5]),
             numpy.random.uniform(size=[3, 4, 5]),
             numpy.random.uniform(size=[3, 6, 5])]

    ilist = [numpy.array([1, 11]),
              numpy.array([2, 3, 4, 8]),
              numpy.array([0, 5, 6, 7, 9, 10]) ]

    out = concatenate_reorder(alist, ilist, axis=1)

    acat = numpy.concatenate(alist, axis=1)
    icat = numpy.concatenate(ilist, axis=0)

    # Inverse the icap permutation
    inv = numpy.copy(icat)
    for i in range(len(inv)):
        inv[i] = numpy.where(icat == i)[0]

    assert (out == acat[:, inv, :]).all()


def test_balanced_cpu_list():
    cpu_cost_list = [[0, 3], [46, 1], [13, 2], [63, 5], [96, 4]]
    node_number = 3
    asserted_result = [[4], [5], [1, 2, 3]]
    result = balanced_cpu_list(cpu_cost_list, node_number)
    print("asserted_result", asserted_result)
    print("result", result)
    assert(asserted_result == result)
