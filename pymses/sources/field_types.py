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
:mod:`pymses.sources.field_types` --- Data source field types
-------------------------------------------------------------

"""
import numpy


class Field(object):
    """
    Data file abstract Field class

    Parameters
    ----------
    name: ``string``
        field name
    ivars: ``list`` of ``int``
        field variable indices in the data file
    dtype: ``numpy.dtype``
        data type.
    """
    _type_name = None

    def __init__(self, name, ivars, dtype):
        self._name = name
        # Ensure integer ivars
        self._ivars = [int(ivar) for ivar in ivars]
        self._dtype = dtype

    def collect_data(self, data):
        """
        Collect field values out of a data array

        Parameters
        ----------
        data: ``numpy.ndarray``

        Returns
        -------
        v: ``numpy.ndarray``
            field values
        """
        return data[..., self._ivars]

    @property
    def variable_indices(self):
        """
        Field variable index list
        """
        return self._ivars

    @property
    def field_name(self):
        """
        Field name
        """
        return self._name

    @property
    def data_type(self):
        return self._dtype

    @classmethod
    def field_type(cls):
        """
        Field type
        """
        return cls._type_name

    def gather(self, data, dset):
        """
        Gather field data in a dataset (abstract method)
        """
        raise NotImplementedError()

    def copy(self):
        """
        Returns
        -------
        f: :class:`~pymses.sources.field_types.Field`
            copy of the itself
        """
        raise NotImplementedError()

    def update_ivars(self, new_ivars):
        """
        Parameters
        ----------
        new_ivars: ``list``of ``int``
            new field variable index list
        """
        if len(new_ivars) != len(self._ivars):
            raise AttributeError("New 'ivars' array must have the same size as the old one.")
        self._ivars = [int(ivar) for ivar in new_ivars]

    def __repr__(self):
        return "%s(name = \"%s\")" % (self.__class__.__name__, self._name)


class Scalar(Field):
    """
    Generic scalar field class describing a 1-dimensional field

    Parameters
    ----------
    name: ``string``
        scalar field name
    ivar: ``int``
        scalar field variable index in the data file
    dtype: ``numpy.dtype``
        data type.
    """
    _type_name = "scalar_field"

    def __init__(self, name, ivar, dtype):
        super(Scalar, self).__init__(name, [ivar], dtype)

    def collect_data(self, data):
        """
        Collect scalar field values out of a data array

        Parameters
        ----------
        data: ``numpy.ndarray``

        Returns
        -------
        v: ``numpy.ndarray``
            scalar field values
        """
        return data[..., self._ivars[0]]

    def copy(self):
        """
        Returns
        -------
        f: :class:`~pymses.sources.field_types.Scalar`
            copy of the Scalar field.
        """
        return Scalar(self._name, self._ivars[0], self._dtype)

    def gather(self, data, dset):
        """
        Gather scalar field data in a dataset

        Parameters
        ----------
        data: ``numpy.ndarray``
        dset: :class:`~pymses.core.datasets.AbstractDataset`
            dataset object
        """
        dset.add_scalars(self._name, self.collect_data(data))

    def __repr__(self):
        return "%s, ivar = %d)" % (super(Scalar, self).__repr__()[:-1], self._ivars[0])


class Vector(Field):
    """
    Generic vector field class describing a n-dimensional field

    Parameters
    ----------
    name: ``string``
        vector field name
    ivars: ``list`` of ``int``
        vector field variable indices in the data file
    dtype: ``numpy.dtype``
        data type.
    """
    _type_name = "vector_field"

    def __init__(self, name, ivars, dtype):
        super(Vector, self).__init__(name, ivars, dtype)

    def copy(self):
        """
        Returns
        -------
        f: :class:`~pymses.sources.field_types.Vector`
            copy of the Vector field.
        """
        return Vector(self._name, self._ivars, self._dtype)

    def gather(self, data, dset):
        dset.add_vectors(self._name, self.collect_data(data))

    def __repr__(self):
        return "%s, ivars = %s)" % (super(Vector, self).__repr__()[:-1], self._ivars)


class MultiValued(Field):
    """
    Generic multivalued field class describing a multivalued field (contiguous variables in the data file)

    Parameters
    ----------
    name: ``string``
        multivalued field name
    ivar_first: ``int``
        multivalued field first variable index in the data file
    nb_vars: ``int``
        multivalued field total number of variables
    dtype: ``numpy.dtype``
        data type.
    """
    _type_name = "multivalued_field"

    def __init__(self, name, ivar_first, nb_vars, dtype):
        self._ifirst = int(ivar_first)
        self._nb_vars = int(nb_vars)
        super(MultiValued, self).__init__(name, list(range(self._ifirst, self._ifirst + self._nb_vars)), dtype)

    def copy(self):
        """
        Returns
        -------
        f: :class:`~pymses.sources.field_types.MultiValued`
            copy of the MultiValued field.
        """
        return MultiValued(self._name, self._ifirst, self._nb_vars, self._dtype)

    def gather(self, data, dset):
        dset.add_multivalued(self._name, self.collect_data(data))

    def update_ivars(self, new_ivars):
        # Assert contiguous variables indices
        if len(new_ivars) > 1:
            if not (numpy.diff(new_ivars) == 1).all():
                raise AttributeError("Multivalued field variable indices must be contiguous. Got %s" % new_ivars)
        super(MultiValued, self).update_ivars(new_ivars)
        self._ifirst = new_ivars[0]

    def __repr__(self):
        return "%s, ivars = [%d, ..., %d])" % (super(MultiValued, self).__repr__()[:-1],
                                               self._ivars[0], self._ivars[-1])


__all__ = ["Scalar", "Vector", "MultiValued"]
