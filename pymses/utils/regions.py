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
:mod:`pymses.utils.regions` --- Regions module
==============================================

"""
import numpy


class Region(object):
    r"""
    Generic region class

    """
    def contains(self, points):
        r"""
        Parameters
        ----------
        points : ``float array of 3D points coordinates``

        Returns
        -------
        points : ``boolean array``
            True when points coordinates are inside the region
        """
        raise NotImplementedError()

    def get_bounding_box(self):
        raise NotImplementedError()

    def filter_points(self, points):
        mask = self.contains(points)
        return points[mask]

    def get_volume(self):
        raise NotImplementedError()

    def random_points(self, npoints, ensure_exact_count=True):
        r"""
        Generates a set of randomly distrubuted points in the region

        Parameters
        ----------
        npoints            : ``int``
            number of points to generate
        ensure_exact_count : ``boolean`` (default True)
            whether the exact required number of random points are generated or not

        Returns
        -------
        points : ``array``
            ramdom points array

        """
        pt_list = []
        npoints_out = 0

        xmin, xmax = self.get_bounding_box()
        delta = xmax - xmin
        ndim = len(xmin)
        nbatch = max(npoints / 4, 1000)

        shape = (nbatch, ndim)

        if not ensure_exact_count:
            filtered_pts = self.filter_points(
                numpy.random.uniform(size=shape, low=0, high=1) * delta + xmin)
            return filtered_pts

        while npoints_out < npoints:
            filtered_pts = self.filter_points(
                numpy.random.uniform(size=shape, low=0, high=1) * delta + xmin)

            npoints_out += filtered_pts.shape[0]
            pt_list.append(filtered_pts)

        return numpy.concatenate(pt_list, axis=0)[:npoints, :]


class Sphere(Region):
    r"""
    Spherical region class

    Parameters
    ----------
    center  : 3-``tuple`` of ``float``
        sphere center coordinates
    radius  : ``float``
        radius of the sphere

    Examples
    --------

        >>> sph = Sphere((0.5, 0.5, 0.5), 1.0)

    """

    def __init__(self, center, radius):
        """
        Spherical region

        """
        self.center = numpy.array(center).ravel()
        self.radius = radius
        self.ndim = len(self.center)
        assert (self.ndim == 3)

    def contains(self, points):
        r"""
        TODO

        """
        assert points.shape[1] == self.ndim

        dx2 = (points - self.center) ** 2
        r2 = numpy.sum(dx2, axis=1)
        return (r2 <= self.radius ** 2)

    def get_bounding_box(self):
        r"""
        TODO

        """
        return (self.center - self.radius, self.center + self.radius)

    def get_volume(self):
        r"""
        Returns
        -------
        V : ``float``
            volume of the sphere (radius :math:`r`) given by :math:`V = \frac{4}{3} \pi r^{3}`

        """
        return 4.0 / 3.0 * numpy.pi * self.radius ** 3


class SphericalShell(Region):
    r"""
    Spherical shell class

    Parameters
    ----------
    center     : 3-``tuple`` of ``float``
        spherical shell center coordinates
    radius_in  : ``float``
        radius of the innerr sphere
    radius_out : ``float``
        radius of the outer sphere

    Examples
    --------

        >>> sph_shell = SphericalShell((0.5, 0.5, 0.5), 0.5, 0.6)

    """

    def __init__(self, center, radius_in, radius_out):
        """ Spherical shell region
        """
        self.center = numpy.array(center).ravel()
        self.radius_in = radius_in
        self.radius_out = radius_out
        self.ndim = len(self.center)
        assert (self.ndim == 3)

    def contains(self, points):
        r"""
        TODO

        """
        assert points.shape[1] == self.ndim

        dx2 = (points - self.center) ** 2
        r2 = numpy.sum(dx2, axis=1)
        return ((r2 <= self.radius_out ** 2) * (self.radius_in ** 2 <= r2))

    def get_bounding_box(self):
        r"""
        TODO

        """
        return (self.center - self.radius_out, self.center + self.radius_out)

    def get_volume(self):
        r"""
        Returns
        -------
        V : ``float``
            volume of the spherical shell (:math:`r_{in}<r<r_{out}`) given by :math:`V = \frac{4}{3} \pi (r_{out}^{3}-r_{in}^{3})`

        """
        return 4.0 / 3.0 * numpy.pi * (self.radius_out ** 3 - self.radius_in ** 3)


class Box(Region):
    r"""
    Box region class

    Parameters
    ----------
    bounds  : 2-``tuple`` of ``list``
        box region boundary min and max positions as a (min, max) ``tuple`` of coordinate arrays

    Examples
    --------
        >>> min_coords = [0.1, 0.2, 0.25]
        >>> max_coords = [0.9, 0.8, 0.75]
        >>> b = Box((min_coords, max_coords))

    """

    def __init__(self, bounds):
        """ Box region

        """
        self.min_coords, self.max_coords = [numpy.asarray(b) for b in bounds]
        self.ndim = len(self.min_coords)

    def contains(self, pointsOrBox):
        if isinstance(pointsOrBox, Box):
            test = ((pointsOrBox.min_coords >= self.min_coords) * (pointsOrBox.max_coords <= self.max_coords))
            return test.all()
        else:
            assert pointsOrBox.shape[1] == self.ndim
            test = ((pointsOrBox >= self.min_coords) * (pointsOrBox <= self.max_coords))
            return test.all(axis=1)

    def __eq__(self, other):
        if not isinstance(other, Box):
            return False
        if not numpy.allclose(other.min_coords, self.min_coords):
            return False
        return numpy.allclose(other.max_coords, self.max_coords)

    def __ne__(self, other):
        return not self.__eq__(other)

    def get_bounding_box(self):
        r"""
        Returns
        -------
        (min_coords, max_coords) : 2-``tuple`` of ``list``
                    bounding box limit
        """
        return (self.min_coords, self.max_coords)

    def get_volume(self):
        r"""
        Returns
        -------
        V : ``float``
            volume of the box given by :math:`V = \displaystyle\prod_{1 \leq i \leq \text{ndim}} (\text{cmax}_{i}-\text{cmin}_{i})`
        """
        return numpy.prod(self.max_coords - self.min_coords)

    def printout(self):
        r"""
            Print bounding box limit in console

        """
        print("Box min", self.min_coords)
        print("Box max", self.max_coords)


class Cube(Box):
    r"""
    Cubic region class

    Parameters
    ----------
    center  : ``tuple``
        cube center coordinates
    width   : ``float``
        size of the cube

    Examples
    --------

        >>> cu = Cube((0.5, 0.5, 0.5), 1.0)

    """

    def __init__(self, center, width):
        """ Cube region
        """
        bounds = numpy.transpose(numpy.array([[coord - width / 2.0, coord + width / 2.0] for coord in center]))
        super(Cube, self).__init__(bounds)
        self.center = center
        self.width = width

    def get_volume(self):
        r"""
        Returns
        -------
        V : ``float``
            volume of the cube (size :math:`L`) given by :math:`V = L^{\text{ndim}}`

        """
        return self.width ** self.ndim


UnitCube = Cube((0.5, 0.5, 0.5), 1.0)


class Cylinder(Region):
    r"""
    Cylinder region class

    Parameters
    ----------
    center       : 3-``tuple`` of ``float``
        cylinder center coordinates
    axis_vector  : 3-``tuple`` of ``float``
        cylinder axis vector coordinates
    radius       : ``float``
        cylinder radius
    height       :	``float``
        cylinder height

    Examples
    --------

        >>> center = (0.5, 0.5, 0.5)
        >>> axis = (0.1, 0.9, -0.1)
        >>> radius = 0.3
        >>> h = 0.05
        >>> cyl = Cylinder(center, axis, radius, h)

    """

    def __init__(self, center, axis_vector, radius, height):
        self.center = numpy.array(center)
        self.axis_vector = numpy.array(axis_vector)
        self.axis_vector /= numpy.linalg.norm(self.axis_vector, ord=2)
        self.height = float(height)
        self.radius = float(radius)

    def contains(self, points):
        r"""
        TODO

        """
        assert points.shape[1] == 3

        p = points - self.center
        u = numpy.dot(p, self.axis_vector)
        pT = p - numpy.outer(u, self.axis_vector)
        v = numpy.sqrt(numpy.sum(pT ** 2, axis=1))

        test_rad = (v <= self.radius)
        test_hgt = (abs(u) <= self.height / 2)
        return test_rad * test_hgt

    def get_bounding_box(self):
        r"""
        TODO

        """
        rad = numpy.sqrt(self.radius ** 2 + (self.height / 2.) ** 2)
        return (self.center - rad, self.center + rad)

    def get_volume(self):
        r"""
        Returns
        -------
        V : ``float``
            volume of the cylinder (radius :math:`r`, height :math:`h`) given by :math:`V = \pi r^{2} h`

        """
        return numpy.pi * self.radius ** 2 * self.height

__all__ = ["Region", "Sphere", "SphericalShell", "Box", "Cube", "UnitCube", "Cylinder"]
