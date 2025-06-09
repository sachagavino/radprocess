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
:mod:`pymses.sources.ramses.field_descr` --- RAMSES output files field descriptor package
-----------------------------------------------------------------------------------------

"""
import numpy as N
from ..field_types import Scalar, Vector, MultiValued


class RamsesFileFieldDescriptor(object):
    """
    Generic Ramses output file field descriptor

    Parameters
    ----------
    field_list : ``list``
         :class:`~pymses.sources.field_types.Field` instance list
    """

    def __init__(self, field_list=None):
        self._field_list = []
        self._fields_by_name = {}

        if field_list is not None:
            for fd in field_list:
                self._add_field(fd)

    @property
    def field_name_list(self):
        """
        Field name list
        """
        return [fd.field_name for fd in self._field_list]

    def __getitem__(self, item):
        fd = self._fields_by_name.get(item)
        if fd is None:
            raise AttributeError("Unknown field '%s'" % item)
        return fd

    def _add_field(self, field):
        """
        Protected method to add a field into the field list

        Parameters
        ----------

        field: Ramses field object
        :raise AttributeError: occurs if the name of the field is already declared into the field list
        """
        name = field.field_name
        if name in self._fields_by_name:
            raise AttributeError("Field '%s' is already declared. Update it or choose a different"
                                 " name for the new field." % name)

        # Gather fields in a file type-base dictionary and a name-based dictionary
        self._field_list.append(field)
        self._fields_by_name[field.field_name] = field

    def add_scalar(self, *args, **kwargs):
        """
        Add a new scalar field into the field list

        Parameters
        ----------
        args: argument list passed along to the `Scalar` constructor
        kwargs: keyword argument dict passed along to the `Scalar` constructor
        """
        self._add_field(Scalar(*args, **kwargs))

    def add_vector(self, *args, **kwargs):
        """
        Add a new vector field into the field list

        Parameters
        ----------
        args: argument list passed along to the `Vector` constructor
        kwargs: keyword argument dict passed along to the `Vector` constructor
        """
        self._add_field(Vector(*args, **kwargs))

    def add_multivalued(self, *args, **kwargs):
        """
        Add a new multivalued field into the field list

        Parameters
        ----------
        args: argument list passed along to the `MultiValued` constructor
        kwargs: keyword argument dict passed along to the `MultiValued` constructor
        """
        self._add_field(MultiValued(*args, **kwargs))

    def remove_field(self, field_name):
        """
        Removes a field with a given name from the field list

        Parameters
        ----------
        field_name: name of the field to remove
        """
        fd = self.__getitem__(field_name)
        self._field_list.remove(fd)
        del self._fields_by_name[field_name]

    def get_file_prefix(self, field_name):
        """
        Get the file prefix of a given field (described by its field name)

        :param field_name: Name of the required field
        :rtype : string
        :return: The file prefix of the selected field
        """
        fd = self.__getitem__(field_name)
        return fd.file_prefix

    def reset(self):
        """
        Reset the field list

        """
        self._field_list = []
        self._fields_by_name = {}

    def copy(self):
        """
        Returns
        -------
        c: : class: `~pymses.sources.ramses.field_descr.RamsesFileFieldDescriptor`
            deep copy of itself
        """
        return self.__class__([fd.copy() for fd in self._field_list])

    def __repr__(self):
        return "%s : [%s]" % (self.__class__.__name__, ",\n".join([fd.__repr__() for fd in self._field_list]))

    def gather_read_fields(self, read_fields):
        """
        Get the (ivars, field) list for a set of field names to read

        Parameters
        ----------

        read_fields: ``list`` of ``string``
            list of field names to read

        Return
        ------
        l : ``list``
            (ivars, field) list for the required set of fields to read.

        """
        ivars_to_read = []
        new_fields = []
        field_list = []

        for fname in read_fields:
            field = self.__getitem__(fname)
            field_list.append(field)

        for field in field_list:
            ivars_to_read += field.variable_indices  # List concatenation of all ivars list of every field

        # Uniquify (+sort) the ivar list
        unique_ivars, rev_indices = N.unique(ivars_to_read, return_inverse=True)

        # Gather a new field list, each of which has its ivars list reordered according to the unique ivar list
        i = 0
        for field in field_list:
            fd_new = field.copy()
            j = len(field.variable_indices)
            fd_new.update_ivars(rev_indices[i:i + j])
            i += j
            new_fields.append(fd_new)

        return [int(i) for i in unique_ivars], new_fields


class _RamsesFilePrefixMixin(object):
    """
    Ramses file prefix mixin
    """
    def _set_prefix(self, file_prefix):
        self._file_prefix = file_prefix

    @property
    def file_prefix(self):
        """
        Ramses file prefix
        """
        return self._file_prefix


class RamsesAmrScalar(Scalar, _RamsesFilePrefixMixin):
    """
    Ramses scalar field descriptor

    Parameters
    ----------

    file_prefix : ``string``
        Ramses output data file prefix
    name : ``string``
        name of the field
    ivar : ``int``
        data file variable index that needs to be read

    Examples
    --------

        >>> rho_field = RamsesAmrScalar("hydro", "rho", 0)
        >>> P_field = RamsesAmrScalar("hydro", "P", 4)

    See also
    --------
    RamsesOutput.define_amr_scalar_field : :meth:`~pymses.sources.ramses.output.RamsesOutput.define_amr_scalar_field`
        lets you define a scalar field into a :class:`~pymses.sources.ramses.output.RamsesOutput` instance.
    """
    def __init__(self, file_prefix, name, ivar):
        super(RamsesAmrScalar, self).__init__(name, ivar, dtype=N.float64)
        self._set_prefix(file_prefix)

    def copy(self):
        """
        Returns
        -------
        f: :class:`~pymses.sources.ramses.field_descr.RamsesAmrScalar`
            copy of the RamsesAmrScalarField
        """
        return RamsesAmrScalar(self._file_prefix, self._name, self._ivars[0])

    def __repr__(self):
        return "%s, file = \"%s\")" % (super(RamsesAmrScalar, self).__repr__()[:-1], self.file_prefix)


class RamsesAmrVector(Vector, _RamsesFilePrefixMixin):
    """
    Ramses vector field descriptor

    Parameters
    ----------

    file_prefix : ``string``
        Ramses output data file prefix
    name : ``string``
        name of the field
    ivars : ``int``
        data file variable index list that needs to be read

    Examples
    --------

        >>> velocity_field = RamsesAmrVector("hydro", "v", [1, 2, 3])

    See also
    --------
    RamsesOutput.define_amr_vector_field : :meth:`~pymses.sources.ramses.output.RamsesOutput.define_amr_vector_field`
        lets you define a vector field into a :class:`~pymses.sources.ramses.output.RamsesOutput` instance.
    """
    def __init__(self, file_prefix, name, ivars):
        super(RamsesAmrVector, self).__init__(name, ivars, dtype=N.float64)
        self._set_prefix(file_prefix)

    def copy(self):
        """
        Returns
        -------
        f: :class:`~pymses.sources.ramses.field_descr.RamsesAmrVector`
            copy of the RamsesAmrVector
        """
        return RamsesAmrVector(self._file_prefix, self._name, self._ivars[:])

    def __repr__(self):
        return "%s, file = \"%s\")" % (super(RamsesAmrVector, self).__repr__()[:-1], self.file_prefix)


class RamsesAmrMultiValued(MultiValued, _RamsesFilePrefixMixin):
    """
    Ramses multi-valued field descriptor (contiguous variables)

    Parameters
    ----------

    file_prefix : ``string``
        Ramses output data file prefix
    name : ``string``
        name of the field
    ivar_first : ``int``
        data file variable index of the first variable that needs to be read
    nb_vars : ``int``
        number of data file variable that needs to be read

    Examples
    --------

        >>> B_field = RamsesAmrMultiValued("hydro", "B", 4, 6)

    See also
    --------
    RamsesOutput.define_amr_multivalued_field : :meth:`~pymses.sources.ramses.output.RamsesOutput.define_amr_multivalued_field`
        lets you define a multivalued field into a :class:`~pymses.sources.ramses.output.RamsesOutput` instance.
    """
    def __init__(self, file_prefix, name, ivar_first, nb_vars):
        super(RamsesAmrMultiValued, self).__init__(name, ivar_first, nb_vars, dtype=N.float64)
        self._set_prefix(file_prefix)

    def copy(self):
        """
        Returns
        -------
        f: :class:`~pymses.sources.ramses.field_descr.RamsesAmrMultiValued`
            copy of the RamsesAmrMultiValued
        """
        return RamsesAmrMultiValued(self._file_prefix, self._name, self._ivars[0], self._ivars[-1] - self._ivars[0] + 1)

    def __repr__(self):
        return "%s, file = \"%s\")" % (super(RamsesAmrMultiValued, self).__repr__()[:-1], self.file_prefix)


class AmrDataFileFieldDescriptor(RamsesFileFieldDescriptor):
    """
    Ramses AMR data file field descriptor

    Parameters
    ----------

    field_list : ``list``
         :class:`~pymses.sources.field_types.Field` instance list
    """
    def __init__(self, field_list=None):
        self._fields_by_file = {}
        super(AmrDataFileFieldDescriptor, self).__init__(field_list)

    def add_scalar(self, *args, **kwargs):
        """
        Add a new AMR scalar field into the field list

        Parameters
        ----------
        args: argument list passed along to the `RamsesAmrScalar` constructor
        kwargs: keyword argument dict passed along to the `RamsesAmrScalar` constructor
        """
        field = RamsesAmrScalar(*args, **kwargs)
        self._add_field(field)

    def add_vector(self, *args, **kwargs):
        """
        Add a new AMR vector field into the field list

        Parameters
        ----------
        args: argument list passed along to the `RamsesAmrVector` constructor
        kwargs: keyword argument dict passed along to the `RamsesAmrVector` constructor
        """
        self._add_field(RamsesAmrVector(*args, **kwargs))

    def add_multivalued(self, *args, **kwargs):
        """
        Add a new AMR multivalued field into the field list

        Parameters
        ----------
        args: argument list passed along to the `RamsesAmrMultiValued` constructor
        kwargs: keyword argument dict passed along to the `RamsesAmrMultiValued` constructor
        """
        self._add_field(RamsesAmrMultiValued(*args, **kwargs))

    def _add_field(self, field):
        """
        Protected method to add a field into the field list

        Parameters
        ----------

        field: Ramses field object
        """
        super(AmrDataFileFieldDescriptor, self)._add_field(field)

        prefix = field.file_prefix
        if prefix not in self._fields_by_file:  # Init file type field field list
            self._fields_by_file[prefix] = [field]
        else:  # Append field to the file type field list
            self._fields_by_file[prefix].append(field)

    def remove_field(self, field_name):
        """
        Removes a field with a given name from the field list and the field-by-file dict.

        Parameters
        ----------
        field_name: name of the field to remove
        """
        fd = self.__getitem__(field_name)
        super(AmrDataFileFieldDescriptor, self).remove_field(field_name)
        flist = self._fields_by_file[fd.file_prefix]
        flist.remove(fd)
        if len(flist) == 0:  # Delete file type field list if empty
            del self._fields_by_file[fd.file_prefix]

    def reset(self):
        """
        Reset the field list

        """
        super(AmrDataFileFieldDescriptor, self).reset()
        self._fields_by_file = {}

    def iter_fields_by_file(self):
        for file_prefix, field_list in self._fields_by_file.items():
            yield file_prefix, field_list

    def gather_read_fields(self, amr_read_fields):
        """
        Get the (ivars, field) amr file type-based dictionary for a set of amr field names to read

        Parameters
        ----------

        amr_read_fields: ``list`` of ``string``
            list of amr field names to read

        Return
        ------
        d : ``dict``
            (ivars, field) amr file type-based dictionary for the required set of amr fields containing 2-tuples of the
            list if the variable indices to read and the list of field.

        """
        descr = AmrDataFileFieldDescriptor([self.__getitem__(fname) for fname in amr_read_fields])

        ivars_descrs_by_file = {}

        for prefix, fds in descr.iter_fields_by_file():
            fname_list = [f.field_name for f in fds]
            ivars, new_fields = super(AmrDataFileFieldDescriptor, descr).gather_read_fields(fname_list)

            # Set the (unique ivar array, new field list) tuple for the current Ramses file prefix
            ivars_descrs_by_file[prefix] = (ivars, new_fields)

        return ivars_descrs_by_file


class SinkParticleFileFieldDescriptor(RamsesFileFieldDescriptor):
    def __init__(self, field_list, position_field_name):
        self._pos_field_name = position_field_name
        super(SinkParticleFileFieldDescriptor, self).__init__(field_list)

    @property
    def sink_position_field(self):
        return self._pos_field_name


__all__ = ["RamsesAmrScalar", "RamsesAmrVector", "RamsesAmrMultiValued", "AmrDataFileFieldDescriptor",
           "SinkParticleFileFieldDescriptor"]
