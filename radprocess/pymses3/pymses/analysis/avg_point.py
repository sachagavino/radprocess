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


def average_point(source, weight_func=None, returned=False):
    r"""
    Return the average point coordinates of a PointDataSource assuming an optional weight function

    Parameters
    ----------

    source : :ref:`PointDataSource`
        the `PointDataSource` from which the average point is computed
    weight_func : ``function``, optional
        `function` used to give a weight for each point of the PointDataSource. Takes a :ref:`AbstractDataset` for single
        argument and returns the weight value for each point
    returned : boolean, optional (default False)
        if True, the sum of the weights is also returned

    Returns
    -------
    av_pos : ``array``
        coordinates of the barycenter
    sow : ``float``
        returned only if `returned` was True. Sum of the weights

    """
    p = []
    sow = []

    for dset in source.iter_dsets():
        if dset.npoints == 0:
            continue
        p_i, sow_i = dset.average_point(weight_func=weight_func)
        p.append(p_i)
        sow.append(sow_i)

    p = numpy.asarray(p)
    sow = numpy.asarray(sow)

    if returned:
        p_tot, sow_tot = numpy.average(p, axis=0, weights=sow, returned=True)
        return (p_tot, sow_tot)
    else:
        p_tot = numpy.average(p, axis=0, weights=sow)
        return p_tot

__all__ = ["average_point"]
