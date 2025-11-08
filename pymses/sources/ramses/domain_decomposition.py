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
:mod:`pymses.sources.ramses.domain_decomposition` --- Base RAMSES domain decomposition package
----------------------------------------------------------------------------------------------

"""
import numpy
from . import hilbert as H


class BaseDomainDecomposition(object):
    def __init__(self):
        super(BaseDomainDecomposition, self).__init__()

    def map_points(self, points):
        """Returns the domain ids for each of the input points.

        Parameters
        ----------
        points : (npoints, ndim) numpy ``array``

        Returns
        -------
        point_domain_list: ``list`` of point mapped domains
        """
        raise NotImplementedError()

    def map_region(self, region):
        """Returns a list of domain ids which ensure covering of the given region

        Parameters
        ----------
            region -- a Region object

        Returns
        -------
        point_domain_list: ``list`` of region mapped domains
        """
        return self.map_box(region.get_bounding_box())

    def map_box(self, box, read_lmax=None):
        """
        Returns a list of all the domains ids which fully cover the given box

        Parameters
        ----------
        box -- ``tuple``
            (min_coords, max_coords) tuple
        read_lmax: maximum read AMR level

        Returns
        -------
        point_domain_list: ``list`` of box mapped domains
        """
        raise NotImplementedError()


class HilbertDomainDecomp(BaseDomainDecomposition):
    """
    Peano-Hilbert decomposition of the cube [0, 1[^ndim
    """

    def __init__(self, ndim, keys_min, keys_max, level_bounds):
        super(HilbertDomainDecomp, self).__init__()
        keys_min = numpy.asarray(keys_min)
        keys_max = numpy.asarray(keys_max)

        assert keys_min[0] == 0.0
        assert keys_min.shape == keys_max.shape

        self.keys_min = keys_min
        self.keys_max = keys_max
        self.ncpu = len(keys_min)
        self.levelmin, self.level_max = level_bounds
        self.ndim = ndim
        self.minimal_grid_list = {}
        self.minimal_order_list = {}
        self.minimal_obounds_list = {}
        self.minimal_overlap_grid_list = {}
        self.minimal_overlap_order_list = {}
        self.minimal_overlap_obounds_list = {}

        self.compute_minimal_domain_descr()

    def compute_minimal_domain_descr(self):
        self.minimal_grid_list = {}
        self.minimal_order_list = {}
        self.minimal_obounds_list = {}
        self.minimal_overlap_grid_list = {}
        self.minimal_overlap_order_list = {}
        self.minimal_overlap_obounds_list = {}

        for idomain in range(self.ncpu):
            hkey_min = self.keys_min[idomain]
            hkey_max = self.keys_max[idomain]

            hilbert_order = 0
            iblocks = numpy.zeros((1,self.ndim), dtype='i')

            # Handle domain boundary hilbert key precision error
            if self.level_max > H._LEVEL_HILBERT_KEY_ERROR:
                hkey_error = (1 << (self.ndim*(self.level_max+1 - H._LEVEL_HILBERT_KEY_ERROR)))
                hkey_min = hkey_min - hkey_error
                hkey_max = hkey_max + hkey_error
                if hkey_min < 0:
                    hkey_min = 0
                last_key = (1 << (self.ndim*(self.level_max+1)))
                if hkey_max > last_key:
                    hkey_max = last_key

            block_list, order_list, overlap_block_list, overlap_order_list = \
                    H.hilbert_overlap((hkey_min, hkey_max), iblocks, hilbert_order, self.level_max, self.ndim)
            block_list = numpy.asarray(block_list, dtype='i')
            order_list = numpy.asarray(order_list, dtype='i')
            overlap_block_list = numpy.asarray(overlap_block_list, dtype='i')
            overlap_order_list = numpy.asarray(overlap_order_list, dtype='i')

            nlev = self.level_max-self.levelmin+1
            # Order sorting
            so = numpy.argsort(order_list)
            order_list = order_list[so]
            order_bounds = numpy.searchsorted(order_list, numpy.arange(self.levelmin, self.level_max+2))
            block_list = H.indices_to_positions(block_list[so,:], order_list)

            so = numpy.argsort(overlap_order_list)
            overlap_order_list = overlap_order_list[so]
            overlap_order_bounds = numpy.searchsorted(overlap_order_list, numpy.arange(self.levelmin, self.level_max+2))
            if overlap_block_list.size > 0:
                overlap_block_list = H.indices_to_positions(overlap_block_list[so,:], overlap_order_list)
            self.minimal_grid_list[idomain] = block_list
            self.minimal_order_list[idomain] = order_list
            self.minimal_obounds_list[idomain] = order_bounds

            self.minimal_overlap_grid_list[idomain] = overlap_block_list
            self.minimal_overlap_order_list[idomain] = overlap_order_list
            self.minimal_overlap_obounds_list[idomain] = overlap_order_bounds
            #print "CPU #%5i : (blocks = %4i, overlap = %6i)"%(idomain+1, len(order_list), len(overlap_order_list))

    def map_points(self, points):
        """Returns the domain ids for each of the input points.

        Parameters
        ----------

        points : (npoints, ndim) numpy ``array``

        Returns
        -------

        """

        # Convert points to indices
        indices = H.positions_to_indices(points, self.level_max+1)
        # Evaluate Hilbert keys
        point_keys = H.compute_hilbert_key(indices, self.level_max+1)

        return numpy.digitize(point_keys, self.keys_min)

    def map_box(self, box, read_lmax=H._LEVEL_HILBERT_KEY_ERROR):
        """
        Returns a list of all the domains ids which fully cover the given box
        """

        pmin, pmax = [numpy.asarray(elem) for elem in box]
        # Some sanity checks
        assert (pmin <= pmax).all()

        cpu_set = set()
        for icpu in range(1, self.ncpu+1):
            idomain = icpu-1
            (gl, ol, noverlap) = self.minimal_domain(idomain, read_lmax)
            # For each cpu we loop over the minimal block decomposition to see if there is an overlap or not
            for block,block_order in zip(gl, ol):
                block_size = 1./2**(block_order+1)
                block_min=numpy.zeros(self.ndim)
                block_max=numpy.zeros(self.ndim)
                for i in range(self.ndim):
                    block_min[i]=block[i]-block_size
                    block_max[i]=block[i]+block_size
                if ((pmin < block_max).all() and (pmax > block_min).all()):
                    # There's definitely an overlap: we need to take into account this cpu
                    # print "icpu",icpu, "pmin", pmin, "pmax", pmax,"block_min",block_min,"block_max",block_max
                    cpu_set.add(icpu)
                    break
        #print "sorted(list(cpu_set)) ",sorted(list(cpu_set))
        return sorted(list(cpu_set))

    def minimal_domain(self, idomain, read_lmax=H._LEVEL_HILBERT_KEY_ERROR):
        "Returns the minimal list of grids of the domain for a given cpu domain"
        # Some sanity checks
        assert idomain in range(self.ncpu)
        if read_lmax is None:
            read_ilevel = self.level_max-self.levelmin
        else:
            read_ilevel = max(min(read_lmax,self.level_max), self.levelmin)-self.levelmin
        if read_ilevel > H._LEVEL_HILBERT_KEY_ERROR:
            read_ilevel = H._LEVEL_HILBERT_KEY_ERROR

        i = self.minimal_obounds_list[idomain][read_ilevel+1]
        ol = self.minimal_order_list[idomain][:i]
        gl = self.minimal_grid_list[idomain][:i,:]

        j = self.minimal_overlap_obounds_list[idomain][read_ilevel]
        k = self.minimal_overlap_obounds_list[idomain][read_ilevel+1]
        ool = self.minimal_overlap_order_list[idomain][j:k]
        ogl = self.minimal_overlap_grid_list[idomain][j:k, :]

        noverlap = ool.size
        if noverlap != 0:
            ol = numpy.concatenate([ol, ool], axis=0)
            gl = numpy.concatenate([gl, ogl], axis=0)
        return gl, ol, noverlap


class PlanarDomainDecomposition(BaseDomainDecomposition):
    def __init__(self, ndim, bound_keys):
        super(PlanarDomainDecomposition, self).__init__()
        self.keys = numpy.asarray(bound_keys)

        if self.keys[0] != 0.0:
            raise AttributeError("Minimum domain decomposition bound key value is not 0.0 !")

        # Decomposition keys are in length unit (ranging between 0.0 and boxlen) => normalising to [0.0, 1.0]
        self.keys /= self.keys[-1]
        self.ncpu = len(self.keys)
        self.ndim = ndim

    def map_points(self, points):
        """
        Maps the input points according to their layer, depending on their last dimension coordinates.
        """
        layer_coord = points[:, self.ndim-1]
        return numpy.clip(numpy.searchsorted(self.keys, layer_coord, side='left'), a_min=1, a_max=self.ncpu)

    def map_box(self, box, read_lmax=None):
        """
        Returns a list of all the layer domain ids which fully cover the given box, base on the box last dimension
        coordinates.

        Parameters
        ----------
        box -- ``tuple``
            (min_coords, max_coords) tuple
        read_lmax: maximum read AMR level

        Returns
        -------
        point_domain_list: ``list`` of box mapped domains
        """
        pmin, pmax = [numpy.asarray(elem) for elem in box]
        coord_min = pmin[-1]
        coord_max = pmax[-1]
        cpu_bounds = numpy.clip(numpy.searchsorted(self.keys, [coord_min, coord_max],
                                                   side='left'), a_min=1, a_max=self.ncpu)
        return list(range(cpu_bounds[0], cpu_bounds[1]+1))


__all__ = ["BaseDomainDecomposition", "PlanarDomainDecomposition", "HilbertDomainDecomp"]
