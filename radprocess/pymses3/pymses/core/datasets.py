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
:mod:`pymses.core.datasets` --- PyMSES generic dataset module
=============================================================

"""
import numpy

from pymses.core.misc import random_sinc_dist
from pymses.utils.misc import concatenate_reorder
from .sources import DataSource
from pymses.utils.hdf5 import HDF5Serializable


class AbstractDataset(DataSource, HDF5Serializable):
    r"""
    Base class for all dataset objects

    """
    _pymses_h5_version = 1

    def source_type(self):
        pass

    def __init__(self):
        self._all_data = {}
        self.vectors = []
        self.scalars = []
        self.multivalued = []
        super(AbstractDataset, self).__init__(data_list=[0])

    @property
    def fields(self):
        r"""
        Dictionary of the fields in the dataset

        """
        return self._all_data

    def __getitem__(self, item):
        """
        Quick access to dataset fields
        """
        return self._all_data[item]

    def _add_field(self, name, data):
        if name in self._all_data:
            raise AttributeError("Cannot have duplicate data name in dataset")
        self._all_data[name] = data

    def add_scalars(self, name, data):
        r"""
        Scalar field addition method

        Parameters
        ----------
        name : ``string``
            human-readable name of the scalar field to add
        data : ``numpy.ndarray``
            raw numpy data array of the new scalar field

        """
        self._add_field(name, data)
        self.scalars.append(name)

    def add_vectors(self, name, data):
        r"""
        Vector field addition method

        Parameters
        ----------
        name : ``string``
            human-readable name of the vector field to add
        data : ``numpy.ndarray``
            raw numpy data array of the new vector field

        """
        self._add_field(name, data)
        self.vectors.append(name)

    def add_multivalued(self, name, data):
        r"""
        Multi-valued field addition method

        Parameters
        ----------
        name : ``string``
            human-readable name of the multi-valued field to add
        data : ``numpy.ndarray``
            raw numpy data array of the new multi-valued field

        """
        self._add_field(name, data)
        self.multivalued.append(name)

    def get_domain_dset(self, idomain, verbose=None):
        """
        Get the dataset... returns itself
        """
        return self

    def iter_dsets(self, verbose=None):
        r"""
        Returns an iterator over itself

        Parameters
        ----------
        verbose: ``bool``
            verbosity boolean flag. Default None.
        """
        return iter([self])

    def _h5_serialize(self, h5group, **kwargs):
        """
        Serialize the abstract dataset into a HDF5 file.

        Parameters
        ----------
        h5group: ``h5py.Group``
            HDF5 group object where to store the abstract dataset object.
        """
        def process(kind, kind_vars):
            group = h5group.create_group(kind)
            for name in kind_vars:
                group.create_dataset(name, data=self[name])

        # Write the data
        process("scalars", self.scalars)
        process("vectors", self.vectors)
        process("multivalued", self.multivalued)

    def _read_fields_from_HDF5(self, h5group):
        dtype_list = [("scalars", self.add_scalars),
                      ("vectors", self.add_vectors),
                      ("multivalued", self.add_multivalued)]

        for kind, add_func in dtype_list:
            group = h5group[kind]
            for name, data in group.items():
                add_func(name, data[...])

    @classmethod
    def _h5_deserialize(cls, h5group, version):
        """
        Read an abstract dataset from a HDF5 file.

        Parameters
        ----------
        h5group: ``h5py.Group``
            Main group to deserialize the AbstractDataset from.
        version: ``int``
            Version of the AbstractDataset class to deserialize.

        Returns
        -------
        dset: :class:`~pymses.core.datasets.AbstractDataset`
            new abstract dataset instance
        """
        # Check AbstractDataset object version written in HDF5 file
        if version == cls._pymses_h5_version:  # Current AbstractDataset version (h5py syntax)
            dset = cls()
            dset._read_fields_from_HDF5(h5group)
        else:
            raise AttributeError("Unknown AbstractDataset version (%d)" % int(version))

        return dset


class PointDataset(AbstractDataset):
    r"""
    Point-based dataset base class

    Parameters
    ----------
    points: list of points
    """
    _pymses_h5_version = 1

    def __init__(self, points):
        super(PointDataset, self).__init__()
        self.points = numpy.asarray(points)
        self.npoints = self.points.shape[0]

    def source_type(self):
        return DataSource.PARTICLE_SOURCE

    def transform(self, xform):
        r"""
        Transform the dataset according to the given `xform` :class:`Transformation<pymses.core.transformations.Transformation>`

        Parameters
        ----------
        xform : :class:`Transformation<pymses.core.transformations.Transformation>`
        """
        # Do nothing to the scalars
        pass

        # Transform the vectors at the points self.points
        for vname in self.vectors:
            self._all_data[vname] = xform.transform_vectors(self._all_data[vname], self.points)

        # Do nothing to the multivalued fields
        pass

        # Transform the points
        self.points = xform.transform_points(self.points)

    @classmethod
    def concatenate(cls, dsets, reorder_indices=None):
        r"""
        Datasets concatenation class method. Return a new dataset

        Parameters
        ----------
        dsets : ``list`` of ``PointDataset``
            list of all datasets to concatenate
        reorder_indices : ``array`` of ``int`` (default to None)
            particles reordering indices

        Returns
        -------
        dset : the new created concatenated ``PointDataset``

        """
        # Prefilter the dsets list to remove any None
        dsets = [ds for ds in dsets if ds is not None]

        if dsets == []:
            return None

        # Get the concatenated points and build the new dataset
        dset_points = [ds.points for ds in dsets]
        new_points = concatenate_reorder(dset_points, reorder_indices, axis=0)
        new_dset = cls(new_points)

        dset0 = dsets[0]

        def process(kind_vars, add_func):
            for name in kind_vars:
                # Get the same data across all datasets into a list
                data_list = [ds[name] for ds in dsets]

                # Concatenate, possibly with reordering
                new_data = concatenate_reorder(data_list, reorder_indices, axis=0)

                # Add to dataset with the proper kind
                add_func(name, new_data)

        # Process the data
        process(dset0.vectors, new_dset.add_vectors)
        process(dset0.scalars, new_dset.add_scalars)
        process(dset0.multivalued, new_dset.add_multivalued)

        return new_dset

    @classmethod
    def add_random_shift(cls, points, size):
        r"""
        Add a random shift to point positions in order to avoid grid alignment
        effect on FFT-convolved (spallted) processed images. The "size" (from CellsToPoints Filter
        and IsotropicExtPointDataset) is needed to know the shift amplitude.
        """
        shift = random_sinc_dist(points.shape, power=2)
        return points + shift * size[:, numpy.newaxis]

    def copy(self):
        c = self.__class__(self.points.copy())
        for key in self.scalars:
            c.add_scalars(key, self[key].copy())
        for key in self.vectors:
            c.add_vectors(key, self[key].copy())
        for key in self.multivalued:
            c.add_multivalued(key, self[key].copy())
        return c

    def reorder_points(self, reorder_indices):
        r"""
        Datasets reorder method. Return a new dataset

        Parameters
        ----------
        reorder_indices : ``array`` of ``int``
            points order indices

        Returns
        -------
        dset : the new created reordered ``PointDataset``

        """
        pts = self.points.copy()
        new_dset = self.__class__(pts[reorder_indices])

        def process(kind_vars, add_func):
            for name in kind_vars:
                add_func(name, self[name][reorder_indices])

        # Process the data
        process(self.vectors, new_dset.add_vectors)
        process(self.scalars, new_dset.add_scalars)
        process(self.multivalued, new_dset.add_multivalued)

        return new_dset

    def filtered_by_mask(self, mask_array):
        r"""
        Datasets filter method. Return a new dataset

        Parameters
        ----------
        mask_array : ``numpy.array`` of ``numpy.bool``
            filter mask

        Returns
        -------
        dset : the new created filtered ``PointDataset``

        """
        if self.npoints == 0:
            # Don't attempt to filter an empty data set
            return self

        new_dset = self.__class__(self.points[mask_array])

        def process(kind_vars, add_func):
            for name in kind_vars:
                data = self[name]
                add_func(name, data[mask_array])

        # Reconstruct data fields with filtered data
        process(self.scalars, new_dset.add_scalars)
        process(self.vectors, new_dset.add_vectors)
        process(self.multivalued, new_dset.add_multivalued)

        return new_dset

    def average_point(self, weight_func=None):
        w = None
        if weight_func is not None:
            w = weight_func(self)

        try:
            p, sow = numpy.average(self.points, axis=0, weights=w, returned=True)
        except:
            p = numpy.zeros(self.points.shape[1])
            sow = 0.

        return p, sow

    def _h5_serialize(self, h5group, **kwargs):
        """
        Serialize the point dataset into a HDF5 file.

        Parameters
        ----------
        h5fg: ``h5py.Group``
            HDF5 group object where to store the point dataset object.
        """
        # Write the points
        h5group.create_dataset("points", data=self.points)

        super(PointDataset, self)._h5_serialize(h5group, **kwargs)

    @classmethod
    def _h5_deserialize(cls, h5group, version):
        """
        Read a point dataset from a HDF5 file.

        Parameters
        ----------
        h5group: ``h5py.Group``
            Main group to deserialize the PointDataset from.
        version: ``int``
            Version of the PointDataset class to deserialize.

        Returns
        -------
        dset: :class:`~pymses.core.datasets.PointDataset`
            new point dataset instance
        """
        # Check PointDataset object version written in HDF5 file
        if version == cls._pymses_h5_version:  # Current PointDataset version (h5py syntax)
            dset = cls(h5group["points"][...])
            dset._read_fields_from_HDF5(h5group)
        else:
            raise AttributeError("Unknown PointDataset version (%d)" % int(version))

        return dset


class IsotropicExtPointDataset(PointDataset):
    r"""
    Extended point dataset class
    """
    def __init__(self, points, sizes=None):
        super(IsotropicExtPointDataset, self).__init__(points)
        if sizes is not None:
            self.add_scalars("size", sizes)

    def get_sizes(self):
        r"""
        Returns
        -------
        sizes : ``array``
            point sizes array

        """
        return self._all_data["size"]


__all__ = ["AbstractDataset", "PointDataset", "IsotropicExtPointDataset"]
