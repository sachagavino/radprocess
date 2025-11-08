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
""" point_utils.pyx -- utility functions to generate arrays of points
"""

import numpy
cimport numpy
cimport cython

cdef extern from "adaptive_gaussian_blur.c":
    void adaptive_gaussian_blur_C(double * map_filtered, double * original_map, int * same_value_pixel_size_map, int map_max_size_i, int map_max_size_j)

cdef extern from "adaptive_gaussian_blur.c":
    void compute_same_value_pixel_size_map_C(int * result_map, double * original_map, int map_max_size_i, int map_max_size_j)

@cython.boundscheck(False)
def corner_points(numpy.ndarray[double, ndim=2] points, numpy.ndarray[double, ndim=1] size):
    cdef char npoints = points.shape[0]
    cdef char ndim = points.shape[1]
    cdef char twotondim = numpy.power(2, ndim)
    cdef char ishift, idim, ipoint
    cdef numpy.ndarray[double, ndim=2] corner_points = numpy.repeat(points, twotondim, axis=0)
    cdef numpy.ndarray[double, ndim=1] xc = numpy.zeros(ndim)

    for ipoint in range(npoints):
        for ind in range(twotondim):
            for idim in range(ndim):
                ishift= (ind >> idim) & 1
                xc[idim] = ishift - 0.5
                corner_points[ipoint+ind, idim] = \
                        corner_points[ipoint+ind, idim] + xc[idim] * size[ipoint]
    return corner_points


def adaptive_gaussian_blur(map, numpy.ndarray[int, ndim=2] same_value_pixel_size_map):
    """
    Function that apply an adaptive gaussian blur given a gaussian filter size map, to improve quality of AMR projected map

    Parameters
    ----------
    map      : ``numpy array of float or double``
        the map to filter
    same_value_pixel_size_map     : ``numpy array of int``
        Same shape as map. For each pixel of map, this array give the size in pixel of the gaussian filter to apply.
        Use compute_same_value_pixel_size_map(map) to get this argument on an AMR partly pixelized map
    Returns
    -------
    result_map : ``numpy array of double``
        map filtered !

    """
    cdef numpy.ndarray[double, ndim=2] _map = numpy.array(map, 'f8')
    cdef int map_max_size_i = map.shape[0]
    cdef int map_max_size_j = map.shape[1]
    cdef numpy.ndarray[double, ndim=2] result_map = numpy.zeros((map_max_size_i, map_max_size_j))
    adaptive_gaussian_blur_C(<double *>result_map.data, <double *>_map.data , <int *>same_value_pixel_size_map.data, <int> map_max_size_i, <int> map_max_size_j)
    return result_map

def compute_same_value_pixel_size_map(map):
    """
    Function that compute a real pixel size (i.e. pixels with same values) on a map, for each
    pixel of the map (aimed to be used on AMR ray tracing map with adaptive_gaussian_blur, to find
    the size of the local gaussian filtering to apply)

    Parameters
    ----------
    map      : ``numpy array of float or double``
        the map to filter
    Returns
    -------
    same_value_pixel_size_map : ``numpy array of int``
        same_value_pixel_size_map (the size is in pixel extension)

    """
    cdef numpy.ndarray[double, ndim=2] _map = numpy.array(map, 'f8')
    cdef int map_max_size_i = map.shape[0]
    cdef int map_max_size_j = map.shape[1]
    cdef numpy.ndarray[int, ndim=2] result_map = numpy.zeros((map_max_size_i, map_max_size_j), 'i')
    compute_same_value_pixel_size_map_C(<int *>result_map.data, <double *>_map.data , <int> map_max_size_i, <int> map_max_size_j)
    return result_map

__all__ = ["corner_points", "adaptive_gaussian_blur", "compute_same_value_pixel_size_map"]
