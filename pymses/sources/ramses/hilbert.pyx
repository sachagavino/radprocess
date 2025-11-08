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
r"""
:mod:`pymses.sources.ramses.hilbert` --- Peano-Hilbert domain decomposition tools
----------------------------------------------------------------------------------

"""

import numpy
cimport numpy
cimport cython

# Define integer types used for index arrays
INT_t = int
ctypedef numpy.int_t cINT_t

CHAR_t = numpy.int8
ctypedef numpy.int8_t cCHAR_t

FL_t = numpy.float64
ctypedef numpy.float64_t cFL_t

# Define Hilbert Key max precision to accept minimal domain mapping error
_LEVEL_HILBERT_KEY_ERROR = 6 # Experimaly found value to avoid problems

# Hilbert states
# Those arrays are in Fortran order, so we reshape them accordingly and .copy()
# them so they are contiguous in memory
_HILBERT_3D_STATE = numpy.reshape([1, 2, 3, 2, 4, 5, 3, 5, 0, 1, 3, 2, 7, 6, 4,
    5, 2, 6, 0, 7, 8, 8, 0, 7, 0, 7, 1, 6, 3, 4, 2, 5, 0, 9, 10, 9, 1, 1, 11,
    11, 0, 3, 7, 4, 1, 2, 6, 5, 6, 0, 6, 11, 9, 0, 9, 8, 2, 3, 1, 0, 5, 4, 6, 7,
    11, 11, 0, 7, 5, 9, 0, 7, 4, 3, 5, 2, 7, 0, 6, 1, 4, 4, 8, 8, 0, 6,10, 6, 6,
    5, 1, 2, 7, 4, 0, 3, 5, 7, 5, 3, 1, 1, 11, 11, 4, 7, 3, 0, 5, 6, 2, 1, 6, 1,
    6, 10, 9, 4, 9, 10, 6, 7, 5, 4, 1, 0, 2, 3, 10, 3, 1, 1, 10, 3, 5, 9, 2, 5,
    3, 4, 1, 6, 0, 7, 4, 4, 8, 8, 2, 7, 2, 3, 2, 1, 5, 6, 3, 0, 4, 7, 7, 2, 11,
    2, 7, 5, 8, 5, 4, 5, 7, 6, 3, 2, 0, 1, 10, 3, 2, 6, 10, 3, 4, 4, 6, 1, 7, 0,
    5, 2, 4, 3], (8 ,2, 12), order='F').astype(CHAR_t)

_HILBERT_3D_STATE_r = numpy.reshape([1, 2, 2, 3, 3, 5, 5, 4, 0, 1, 3, 2, 6, 7, 5,
    4, 2, 0, 0, 8, 8, 7, 7, 6, 0, 2, 6, 4, 5, 7, 3, 1, 0, 1, 1, 9, 9, 11, 11, 10,
    0, 4, 5, 1, 3, 7, 6, 2, 11, 6, 6, 0, 0, 9, 9, 8, 3, 2, 0, 1, 5, 4, 6, 7, 9, 7,
    7, 11, 11, 0, 0, 5, 5, 7, 3, 1, 0, 2, 6, 4, 10, 8, 8, 6, 6, 4, 4, 0, 6, 2, 3,
    7, 5, 1, 0, 4, 3, 11, 11, 5, 5, 1, 1, 7, 3, 7, 6, 2, 0, 4, 5, 1, 4, 9, 9, 10,
    10, 6, 6, 1, 5, 4, 6, 7, 3, 2, 0, 1, 5, 10, 10, 1, 1, 3, 3, 9, 6, 4, 0, 2, 3,
    1, 5, 7, 7, 4, 4, 2, 2, 8, 8, 3, 5, 1, 0, 4, 6, 2, 3, 7, 8, 5, 5, 7, 7, 2, 2,
    11, 6, 7, 5, 4, 0, 1, 3, 2, 6, 3, 3, 4, 4, 10, 10,
    2, 3, 1, 5, 7, 6, 4, 0, 2], (8 ,2, 12), order='F').astype(CHAR_t)


_HILBERT_2D_STATE = numpy.reshape([1, 0, 2, 0, 0, 1, 3, 2, 0, 3, 1, 1, 0, 3, 1,
    2, 2, 2, 0, 3, 2, 1, 3, 0, 3, 1, 3, 2, 2, 3, 1, 0], (4, 2, 4),
    order='F').astype(CHAR_t)


# ndim ---> Hilbert state mappings
_HILBERT_NDIM_STATEMAP = {
        2 : _HILBERT_2D_STATE,
        3 : _HILBERT_3D_STATE }
_HILBERT_NDIM_STATEMAP_r = {
        3 : _HILBERT_3D_STATE_r }


def positions_to_indices(pos_array, order):
    """
    Converts an array of positions in [0, 1[^ndim to an array of indices at
    order `order` suitable for compute_hilbert_key
    """
    pos_array = numpy.asarray(pos_array)
    return (pos_array * (1 << order)).astype(INT_t)


def indices_to_positions(ind_array, order):
    """
    Converts an array of indices in [0, 2**order[^ndim to an array of positions in
    [0, 1[^ndim (grid centers)
    """
    ind_array = numpy.asarray(ind_array, 'd')
    return ((ind_array + 0.5) / (1 << order[:,numpy.newaxis]))


@cython.boundscheck(False)
def compute_hilbert_key(ind_array_in, int order):
    """
    Compute the Hilbert keys of points stored in ind_array.

    Arguments:

        ind_array_in -- a (npoints, ndim) numpy array of cell indices.
            See positions_to_indices() to convert floating point positions to an
            index array.

        order -- the order of the Hilbert curve

    Returns:
        a (npoints,) array of integers containing the Hilbert keys
    """
    cdef numpy.ndarray[cINT_t, ndim=2] ind_array = numpy.asarray(ind_array_in, dtype=INT_t)
    cdef int npoints = ind_array.shape[0]
    cdef int ndim = ind_array.shape[1]
    cdef double d_twotondim = float(1 << ndim)

    # Retrieve the state diagram
    cdef numpy.ndarray[cCHAR_t, ndim=3] state_diag
    try:
        state_diag = _HILBERT_NDIM_STATEMAP[ndim]
    except KeyError:
        raise NotImplementedError(
            "ndim=%i unsupported in compute_hilbert_key" % ndim)

    # Bits of the indices for a given point
    cdef numpy.ndarray[cINT_t, ndim=2] ind_bits \
            = numpy.empty([order, ndim], INT_t)

    # Digits of the resulting Hilbert key
    cdef numpy.ndarray[cINT_t, ndim=1] h_digits \
            = numpy.empty(order, INT_t)

    # Hilbert keys in float64
    cdef numpy.ndarray[cFL_t, ndim=1] point_keys \
            = numpy.empty(npoints, FL_t)
    cdef int ipoint, idim, ibit
    cdef int cur_state, new_state, s_digit

    for ipoint in range(npoints):

        # Set ind_bits for current point
        for ibit in range(order):
            for idim in range(ndim):
                ind_bits[ibit, idim] = (ind_array[ipoint, idim] >> ibit) & 1

        # Compute Hilbert key bits
        cur_state = 0
        for ibit in range(order-1, -1, -1):

            # Compute s_digit by interleaving bits
            s_digit = 0
            for idim in range(ndim):
                s_digit += (ind_bits[ibit, idim]) << (ndim-1-idim)

            # Compute the new state from the state diagram
            new_state = state_diag[s_digit, 0, cur_state]
            h_digits[ibit] = state_diag[s_digit, 1, cur_state]
            cur_state = new_state

        # Assemble the point's key as a float 128 bits
        point_keys[ipoint] = 0.0
        for ibit in range(order):
            point_keys[ipoint] += (d_twotondim)**ibit * h_digits[ibit]
    return point_keys


@cython.boundscheck(False)
def compute_indices(hilbert_key_array_in, int order, int ndim):
    """
    Compute the Hilbert keys of points stored in ind_array.

    Arguments:

        hilbert_key_array_in -- a (nkey) numpy array of cell hilbert keys (float128 or uint64).

        order -- the order of the Hilbert curve

        ndim -- dimension of the domain

    Returns:
        a (nkeys, ndim) array of integers containing the cell indices along each dimension
    """
    # Hilbert keys in float64
    cdef numpy.ndarray[cFL_t, ndim=1] point_keys = numpy.asarray(hilbert_key_array_in, FL_t)
    cdef int nkeys = point_keys.size
    cdef double d_twotondim = float(1 << ndim)

    cdef numpy.ndarray[cINT_t, ndim=2] ind_array = \
            numpy.zeros((nkeys, ndim), INT_t)

    # Retrieve the state diagram
    cdef numpy.ndarray[cCHAR_t, ndim=3] state_diag
    try:
        state_diag = _HILBERT_NDIM_STATEMAP_r[ndim]
    except KeyError:
        raise NotImplementedError(
            "ndim=%i unsupported in compute_indices" % ndim)

    # Bits of the indices for a given point
    cdef numpy.ndarray[cINT_t, ndim=2] ind_bits \
            = numpy.empty([order, ndim], INT_t)

    # Digits of the Hilbert keys
    cdef numpy.ndarray[cINT_t, ndim=1] h_digits \
            = numpy.empty(order, INT_t)

    cdef int ikey, idim, ibit
    cdef int cur_state, new_state, s_digit
    cdef cFL_t hkey

    for ikey in range(nkeys):

        # Compute Hilbert key bits
        hkey = point_keys[ikey]
        for ibit in range(order):
            h_digits[ibit] = int((hkey / (d_twotondim)**ibit)%d_twotondim)
            hkey -= h_digits[ibit]*(d_twotondim)**ibit


        # Compute indices bits
        cur_state = 0
        for ibit in range(order-1, -1, -1):

            # Compute the new s_digit from the state diagram
            new_state = state_diag[h_digits[ibit], 0, cur_state]
            s_digit = state_diag[h_digits[ibit], 1, cur_state]
            cur_state = new_state

            # Compute ind_bitd
            for idim in range(ndim):
                ind_bits[ibit, idim] = (s_digit >> (ndim-1-idim)) & 1


        # Set indices for current key
        for ibit in range(order):
            for idim in range(ndim):
                ind_array[ikey, idim] += ind_bits[ibit, idim] << ibit

    return ind_array


def hilbert_overlap(bound_keys, iblocks, hilbert_order, lmax, ndim):
    block_keys = compute_hilbert_key(iblocks, hilbert_order)
    hkey_min, hkey_max = bound_keys
    block_list = []
    order_list = []
    overlap_block_list = []
    overlap_order_list = []

    # 1D resolution at hilbert_order
    nblocks_1d = 2L**hilbert_order

    # This is the ratio of the max Hilbert key in the sim (ie at level
    # level_max+1) to the max Hilbert key at level hilbert_order
    block_key_ratio = (2.0**(lmax+1) / nblocks_1d)**ndim

    # Convert to (level_max+1) min and max Hilbert keys
    fine_keys_min = block_keys * block_key_ratio
    fine_keys_max = (block_keys+1) * block_key_ratio

    # Construct the position tuples of the covered blocks
    gray_blocks = [ numpy.array(
        [ ((bindex >> idim) & 1) for idim in range(ndim)], INT_t )
            for bindex in range(1 << ndim) ]

    for iblock in range(block_keys.size):
        k_min = fine_keys_min[iblock]
        k_max = fine_keys_max[iblock]
        if (hkey_min < k_max) and (hkey_max > k_min):# The block intersects the hilbert domain
            if (hkey_min <= k_min) and (hkey_max >= k_max):
                # Block entirely included in the cpu domain
                block_list.append(iblocks[iblock])
                order_list.append(hilbert_order)
            else:
                # Overlapping block
                overlap_block_list.append(iblocks[iblock])
                overlap_order_list.append(hilbert_order)

                ind = iblocks[iblock]
                niblocks = gray_blocks + (ind[numpy.newaxis, :] << 1)
                if hilbert_order < _LEVEL_HILBERT_KEY_ERROR:
                    bl, ol, ovbl, ovol = hilbert_overlap(bound_keys, niblocks, hilbert_order+1, lmax, ndim)
                    block_list+=bl
                    order_list+=ol
                    overlap_block_list+=ovbl
                    overlap_order_list+=ovol

    return block_list, order_list, overlap_block_list, overlap_order_list


__all__ = ['compute_hilbert_key', 'compute_indices', 'hilbert_overlap', 'indices_to_positions', 'positions_to_indices']
