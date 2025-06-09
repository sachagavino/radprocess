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
:mod:`pymses.analysis.profile_binners` --- Histogram computation module
-----------------------------------------------------------------------

"""

import numpy
from numpy.linalg import norm


class ProfileBinner(object):
    """
    Base class for all profile binners.

    """
    def __init__(self, profile_func, bin_bounds, divide_by_counts=True):
        self.profile_func = profile_func
        self.bin_bounds = bin_bounds
        self.divide_by_counts = divide_by_counts

    def bin_func(self, point_dset):
        """Compute the binned coordinate array for a given PointDataset object
        """
        raise NotImplementedError

    def process(self, source):
        """Compute the profile of the specified data source
        """

        # Prepare full profile histogram
        profile = numpy.zeros(len(self.bin_bounds) - 1)

        for dset in source.iter_dsets():

            bin_coords = self.bin_func(dset)

            # Compute profile for this batch
            dprofile = numpy.histogram(
                bin_coords,
                weights=self.profile_func(dset),
                bins=self.bin_bounds,
                normed=False)[0]

            if self.divide_by_counts:
                # Divide by counts
                counts = numpy.histogram(
                    bin_coords,
                    bins=self.bin_bounds,
                    normed=False)[0]
                counts[counts == 0] = 1
                dprofile = dprofile / counts

            profile += dprofile

        return profile


class SphericalProfileBinner(ProfileBinner):
    """
    Spherical profile binner class

    """

    def __init__(self, center, profile_func, bin_bounds,
                 divide_by_counts=False):
        self.center = numpy.asarray(center)
        super(SphericalProfileBinner, self).__init__(profile_func, bin_bounds, divide_by_counts)

    def bin_func(self, point_dset):
        """Returns the array of distances from `point_dset.points` to
        `self.center` for :class:`PointDataset` objects.
        """
        # Radial vector from center to point_dset.points
        rad = point_dset.points - self.center[numpy.newaxis, :]

        # The bin is determined by the norm of rad
        return numpy.sqrt(numpy.sum(rad * rad, axis=1))


class CylindricalProfileBinner(ProfileBinner):
    """
    Cylindrical profile binner class

    """

    def __init__(self, center, axis_vect, profile_func, bin_bounds, divide_by_counts=False):
        self.center = numpy.asarray(center)
        self.axis_vect = numpy.asarray(axis_vect) / norm(axis_vect, 2)

        super(CylindricalProfileBinner, self).__init__(profile_func, bin_bounds, divide_by_counts)

    def bin_func(self, point_dset):
        """Returns the array of distances from `point_dset.points` to
        the cylinder axis for :class:`PointDataset` objects.
        """
        # Decompose the vector from self.center to point_dset.points into its
        # component along the cylinder axis, and its component in the normal
        # plane
        rad = point_dset.points - self.center[numpy.newaxis, :]
        along = numpy.dot(rad, self.axis_vect)[:, numpy.newaxis] * self.axis_vect
        ortho = rad - along

        # Bin is determined by the radial component only
        return numpy.sqrt(numpy.sum(ortho * ortho, axis=1))


def bin_spherical(source, center, profile_func, bin_bounds, divide_by_counts=False):
    r"""
    Spherical binning function for profile computing

    Parameters
    ----------
    center : ``array``
        center point for the profile
    profile_func : ``function``
        a function taking a ``PointDataset`` object as an input and producing a
        numpy array of weights.
    bin_bounds : ``array``
        a numpy array delimiting the profile bins (see	numpy.histogram documentation)
    divide_by_counts : ``boolean`` (default False)
        if True, the returned profile is the ``array`` containing the sum of weights in each bin.
        if False, the mean weight per bin ``array`` is returned.

    Returns
    -------
    profile : ``array``
        computed spherical profile

    """
    binner = SphericalProfileBinner(center, profile_func, bin_bounds, divide_by_counts)
    return binner.process(source)


def bin_cylindrical(source, center, axis_vect, profile_func, bin_bounds, divide_by_counts=False):
    r"""
    Cylindrical binning function for profile computing

    Parameters
    ----------
    center : ``array``
        center point for the profile
    axis_vect : ``array``
        the cylinder axis coordinates array.
    profile_func : ``function``
        a function taking a ``PointDataset`` object as an input and producing a
        numpy array of weights.
    bin_bounds : ``array``
        a numpy array delimiting the profile bins (see	numpy.histogram documentation)
    divide_by_counts : ``boolean`` (default False)
        if True, the returned profile is the array containing the sum of weights in each bin.
        if False, the mean weight per bin array is returned.

    Returns
    -------
    profile : ``array``
        computed cylindrical profile

    """
    binner = CylindricalProfileBinner(center, axis_vect, profile_func, bin_bounds, divide_by_counts)
    return binner.process(source)


__all__ = ["bin_cylindrical", "bin_spherical"]
