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
:mod:`pymses.core.transformations` Geometrical transformations module
=====================================================================

"""

import numpy
from scipy.linalg import expm
from functools import reduce


class Transformation(object):
    """ Base class for all geometric transformations acting on Numpy arrays
    """

    ndim = 3
    is_affine = False

    # abstract functions ------------------------------------------

    def transform_points(self, coords):
        """Abstract method. Returns transformed coordinates.

        Parameters:
            coords -- a Numpy array with data points along axis 0 and
                coordinates along axis 1+
        """
        raise NotImplementedError()

    def transform_vectors(self, vectors, coords):
        """Abstract method. Returns transformed vector components for vectors
        attached to the provided coordinates.

        Parameters:

            vectors -- a Numpy array of shape (ndata, ndim) containing the
                vector components

            coords -- a Numpy array of shape (ndata, ndim) containing the point
                coordinates
        """
        raise NotImplementedError()

    def inverse(self):
        """Returns the inverse transformation
        """
        raise NotImplementedError()

    # methods -----------------------------------------------------

    def __mul__(self, other):
        """ Overrides multiplication to support intuitive operation composition:
        (T2*T1)(x) = T2(T1(x))
        """
        # T2 is other, T1 is self
        return ChainTransformation([other, self])


class ChainTransformation(Transformation):
    """ Defines the composition of a list of transformations
    """

    def __init__(self, xform_seq):
        """Chains the transformations in the sequence xform_seq to form a new
        transformation.

        NOTE: The transformations are applied in the order of the sequence, ie:
        ChainTransformation([T1, T2]) corresponds to x -> T2(T1(x))
        """

        super(ChainTransformation, self).__init__()

        xform_list = list(xform_seq)
        assert len(xform_list) > 0

        # Check that the dims of the transformations are compatible
        ndim0 = xform_list[0].ndim
        ndims_ok = [xform.ndim == ndim0 for xform in xform_list]
        assert all(ndims_ok)

        # if all transformations are affine, the chained xform is affine
        self.is_affine = all(xform.is_affine for xform in xform_list)
        self.ndim = ndim0
        self.transform_chain = xform_list

    def inverse(self):
        "Inverse of a chained transformation"
        return ChainTransformation(
            [t.inverse() for t in reversed(self.transform_chain)])

    def transform_points(self, coords):
        """
        Applies a chained transformation to coordinates

        :param coords:
        :return:
        """
        # this is the definition of chaining (!)
        return reduce(
            lambda xcoords, xform: xform.transform_points(xcoords),
            self.transform_chain, coords)

    def transform_vectors(self, vectors, coords):
        """
        Applies a chained transformation to vectors

        :param vectors:
        :param coords:
        :return:
        """
        if self.is_affine:
            # Don't transform the coordinates to save time
            xform_coords = lambda pts, xform: pts
        else:
            # We need to transform the coordinates
            xform_coords = lambda pts, xform: xform.transform_points(pts)

        # this is the chain rule for derivatives
        cur_vec = vectors
        cur_pts = coords
        for xform in self.transform_chain:
            cur_vec = xform.transform_vectors(cur_vec, cur_pts)
            cur_pts = xform_coords(cur_pts, xform)

        return cur_vec


class LinearTransformation(Transformation):
    """ A generic (matrix-based) linear transformation
    """
    is_affine = True

    def __init__(self, matrix):
        super(LinearTransformation, self).__init__()

        matrix = numpy.matrix(matrix)
        ndim = matrix.shape[0]
        assert matrix.shape == 2 * (ndim,)

        self.ndim = ndim
        self.matrix = matrix

    def transform_points(self, coords):
        "Applies a linear transformation to coordinates"
        # return the tensor product
        return numpy.tensordot(coords, self.matrix, axes=[1, 1])

    def transform_vectors(self, vectors, coords):
        "Applies a linear transformation to vectors"
        # For a linear transformation, it's the same as transform_points applied
        # to vectors
        return self.transform_points(vectors)

    def inverse(self):
        "Inverse of the linear transformation"
        return LinearTransformation(numpy.linalg.inv(self.matrix))


class AffineTransformation(Transformation):
    """ An affine transformation (of the type x -> L(x) + shift)
    """

    def __init__(self, lin_xform, shift):
        super(AffineTransformation, self).__init__()
        self.shift = numpy.copy(numpy.ravel(shift))
        self.linear_transform = lin_xform
        assert shift.shape == (lin_xform.ndim,)

    def inverse(self):
        "Inverse of an affine transformation"
        inv_lin = self.linear_transform.inverse()
        inv_lin_shift = -inv_lin.transform_points(numpy.array([self.shift])).squeeze()
        return AffineTransformation(inv_lin, inv_lin_shift)

    def transform_points(self, coords):
        "Apply the affine transformation to coordinates"
        return self.linear_transform.transform_points(coords) + self.shift

    def transform_vectors(self, vectors, coords):
        "Apply the affine transformation to vectors"
        return self.linear_transform.transform_vectors(vectors, coords)


# Define the rotation generators
_SO3_GENERATORS = numpy.array([
    [[0, 0, 0],
     [0, 0, -1],
     [0, 1, 0]],
    [[0, 0, 1],
     [0, 0, 0],
     [-1, 0, 0]],
    [[0, -1, 0],
     [1, 0, 0],
     [0, 0, 0]]])


def _ensure_vect(u, ndim=None):
    """ Makes sure that u can be cast into a vector, and returns a sane copy
    of u as a Numpy array if possible

    If ndim if specified, enforces dim(u) == ndim
    """
    uu = numpy.copy(numpy.ravel(u))

    if ndim is not None:
        assert len(uu) == ndim

    return uu


def identity(n):
    """

    Returns
    -------
    the identity as a LinearTransformation object

    """
    m = numpy.eye(n)
    return LinearTransformation(m)


def translation(vect):
    """
    Returns
    -------
    an AffineTransformation object corresponding to a translation
    of the specified vector

    """
    u = numpy.copy(numpy.ravel(vect))
    n = len(u)
    return AffineTransformation(identity(n), u)


def rot3d_axvector_matrix(axis_vect, angle):
    """ Returns the rotation matrix of the rotation with the specified axis
    vector and angle
    """
    # cleanup and normalize the axis vector
    n = _ensure_vect(axis_vect, ndim=3)
    n = n / numpy.linalg.norm(n, ord=2)

    # Construct the rotation matrix from the SO(3) generators
    rotmatrix = expm(angle * numpy.tensordot(_SO3_GENERATORS, n, axes=(0, 0)))

    return numpy.matrix(rotmatrix)


def rot3d_axvector(axis_vect, angle, rot_center=None):
    """ Returns the Transformation corresponding to the rotation specified by
    its axis vector, angle, and rotation center.

    If rot_center is not specified, it is assumed to be [0, 0, 0].
    """
    # first compute the rotation matrix
    rotmatrix = rot3d_axvector_matrix(axis_vect, angle)

    if rot_center is None:
        out = LinearTransformation(rotmatrix)
    else:
        c = _ensure_vect(rot_center, ndim=3)

        # first translate c to 0
        T1 = translation(-c)
        # then rotate
        T2 = LinearTransformation(rotmatrix)
        # then translate back
        T3 = translation(c)

        # Assemble the transformations together in the right order
        out = ChainTransformation([T1, T2, T3])

    return out


def rot3d_euler(axis_sequence, angles, rot_center=None):
    """ Returns the Transformation corresponding to the rotation specified by
    its Euler angles and the corresponding axis sequence convention.

    The rotation is performed by successively rotating the object around its
    current local axis axis_sequence[i] with an angle angle[i], for i = 0, 1, 2.

    See http://en.wikipedia.org/wiki/Euler_angles for details.
    """

    assert len(axis_sequence) == len(angles) == 3

    # initially, local axes coincide with global axes for no transformation
    axloc = numpy.eye(3)
    T = []

    for i in range(3):
        # build the current rotation
        r = rot3d_axvector(axloc[axis_sequence[i], :], angles[i])

        # append it to the list of transformations
        T.append(r)

        # update the current axes in global coords
        axloc = r.transform_points(axloc)

    # add the translation if needed
    if rot_center is not None:
        c = _ensure_vect(rot_center, ndim=3)
        T = [translation(-c)] + T + [translation(c)]

    return ChainTransformation(T)


def rot3d_align_vectors(source_vect, dest_vect, dest_vect_angle=0.0, rot_center=None):
    """
    Gives a :class:`Transformation` which brings a given `source_vect` in alignment
    with a given `dest_vect`.

    Optionally, a second rotation around `dest_vect` can be specified by the parameter `dest_vect_angle`.

    Parameters
    ----------
    source_vect : ``array``
        source vector coordinates array
    dest_vect : ``array``
        destination vector coordinates array
    dest_vect_angle : ``float`` (default 0.0)
        optional final rotation angle around the `dest_vect` vector

    Returns
    -------
    rot : :class:`Transformation`
        rotation bringing `source_vect` in alignment with `dest_vect`.
        This is done by rotating around the normal to the (`source_vect`, `dest_vect`) plane.

    Examples
    --------

        >>> R = rot3d_align_vectors(array([0.,0.,1.]), array([0.5,0.5,0.5]))
    """

    # Normalize input vectors
    source_vect = _ensure_vect(source_vect, 3) / numpy.linalg.norm(source_vect, 2)
    dest_vect = _ensure_vect(dest_vect, 3) / numpy.linalg.norm(dest_vect, 2)

    # Compute the vector to use as the rotation axis
    rot_axis_vect = numpy.cross(source_vect, dest_vect)

    # Compute the angle and normalize the rotation axis
    sin_angle = numpy.linalg.norm(rot_axis_vect, 2)
    cos_angle = numpy.dot(source_vect, dest_vect)
    angle = numpy.arctan2(sin_angle, cos_angle)

    assert sin_angle != 0.0

    rot_axis_vect = rot_axis_vect / sin_angle

    # Create the rotation bringing source_vect onto dest_vect
    rot = rot3d_axvector(rot_axis_vect, angle, rot_center)

    # If needed, do a second rotation around dest_vect
    if dest_vect_angle != 0.0:
        rot = rot3d_axvector(dest_vect, dest_vect_angle) * rot

    return rot


def scale(n, scale_factor, scale_center=None):
    r"""

    """
    m = numpy.eye(n) * scale_factor
    T = LinearTransformation(m)

    if scale_center is not None:
        c = _ensure_vect(scale_center, ndim=n)
        T = ChainTransformation([translation(-c)] + T + [translation(c)])

    return T


__all__ = ["Transformation",
           "AffineTransformation",
           "LinearTransformation",
           "ChainTransformation",
           "identity",
           "translation",
           "rot3d_axvector_matrix",
           "rot3d_axvector",
           "rot3d_euler",
           "rot3d_align_vectors",
           "scale"]
