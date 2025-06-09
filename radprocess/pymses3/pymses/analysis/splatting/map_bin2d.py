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
:mod:`pymses.analysis.splatting.map_bin2d` --- 2D histrogram computation module
-------------------------------------------------------------------------------

"""
import numpy as N


def histo2D(uvw, bins, xyz_range, weights):
    """
    2D binning function

    Parameters
    ----------
    uvw:
    bins:
    xyz_range:
    weights:

    Returns
    -------
    hist:
    """
    # Size/range/pixel size of the histogram
    nx, ny = bins
    xr, yr, zrange = xyz_range
    dx = (xr[1] - xr[0]) / nx
    dy = (yr[1] - yr[0]) / ny

    # Mask inside points
    mask = ((uvw[:, 0] >= xr[0]) * (uvw[:, 0] <= xr[1]) * (uvw[:, 1] >= yr[0]) * (uvw[:, 1] <= yr[1]))
    u = uvw[mask, 0]
    v = uvw[mask, 1]

    # Total number of points to bin
    np = u.size
    if np == 0:
        return dict.fromkeys(weights, N.zeros((nx, ny)))

    # Compute the bin number each point falls into (x an y-axis)
    xcount = N.floor((u - xr[0]) / dx).astype('i')
    ycount = N.floor((v - yr[0]) / dy).astype('i')

    # Using floor, points that fall on an edge are put in the right/upper pixel.
    # For the right/uppermost pixel, we want points on the right/up
    # edge to be counted in the last pixel, and not as a outliers.
    # Rounding precision
    decx = int(-N.log10(dx)) + 6
    decy = int(-N.log10(dy)) + 6
    # Find which points are on the right/uppermost edge
    on_edgex = (N.around(u, decx) == N.around(xr[1], decx))  # [0]
    on_edgey = (N.around(v, decy) == N.around(yr[1], decy))  # [0]
    # Shift these points one bin to the left/bottom.
    xcount[on_edgex] -= 1
    ycount[on_edgey] -= 1

    # 1D index of the points
    xy = N.zeros(np, int)
    xy = xcount * ny + ycount

    binmaps = {}
    for key in weights:
        we = weights[key][mask]

        # Number of repetitions in xy
        flatcount = N.bincount(xy, we)

        # 1D histogram
        hist = N.zeros(nx * ny, float)
        hist[N.arange(len(flatcount))] = flatcount
        binmaps[key] = hist.reshape(nx, ny)

    return binmaps


__all__ = ["histo2D"]
