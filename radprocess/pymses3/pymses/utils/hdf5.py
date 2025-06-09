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
:mod:`pymses.utils.hdf5`: HDF5 serializable interface
-----------------------------------------------------
"""
from .filename import FileUtil

try:
    import h5py
    _h5py_available = True
except ImportError as imperr:
    _h5py_available = False


class HDF5Serializable(object):
    __fact_subclasses_dict = {}
    _pymses_h5_version = 1

    def _h5_serialize(self, h5group, **kwargs):
        raise NotImplementedError()

    def save_HDF5(self, h5fg, **kwargs):
        """
        Returns an serialisable object from a HDF5 file.

        Parameters
        ----------
        h5fg: ``str`` or ``h5py.File`` or ``h5py.Group``
            HDF5 (h5py) file/group or filename.

        Raises
        ------
        imperr: ``ImportError``
            raises an ImportError if the 'h5py' module is not found.
        atterr: ``AttributeError``
            raises an AttributeError if any attribute is not valid.
        ioerr: ``IOError``
            raises an IOError if an error occured while writing into the HDF5 file.
        """
        if not _h5py_available:
            raise ImportError("h5py module is required to handle HDF5 files in PyMSES !")

        if isinstance(h5fg, h5py.Group) and h5fg.__class__ == h5py.Group:
            group = h5fg
            close_when_done = False
            f = None
        else:
            if "where" not in kwargs:
                where = "/"
            else:
                where = kwargs.pop("where")
            f, group, close_when_done = self.__open_h5file(h5fg, where, create=True)

        try:
            # Write version number
            self._h5_version_write(group)

            # Write class name
            self._h5_classname_write(group)

            # Serialize object
            self._h5_serialize(group, **kwargs)
        except Exception as exc:
            raise IOError("HDF5 I/O error : %s" % exc)
        finally:
            if close_when_done and f is not None:
                f.close()

    @classmethod
    def __open_h5file(cls, h5f, where, create=False):
        """
        Returns a HDF5 tables.File object

        Parameters
        ----------
        h5f: ``str`` or ``h5py.File``
            HDF5 (h5py) File object or filename
        where: ``string``
            location of the object in the HDF5 tree
        create: ``bool``
            Open file in writing mode ?

        Returns
        -------
        f: ``h5py.File``
            HDF5 file object.
        g: ``h5py.Group``
            HDF5 group object.

        Raises
        ------
        imperr: ``ImportError``
            raises an ImportError if the 'h5py'  module is not found.
        atterr: ``AttributeError``
            raises an AttributeError if any attribute is not valid.
        """
        # Check the 'where' kwarg is a valid string.
        if not isinstance(where, str):
            raise AttributeError("'where' attribute must be a string. Got '%s'." % type(where))

        if isinstance(h5f, str):
            if create:
                f = h5py.File(FileUtil.new_filepath(h5f, append_extension=FileUtil.HDF5_FILE), 'w')
            else:  # Read an existing HDF5 file
                f = h5py.File(FileUtil.valid_filepath(h5f, append_extension=FileUtil.HDF5_FILE), 'r')

            close_when_done = True
        elif isinstance(h5f, h5py.File):
            f = h5f
            close_when_done = False
        else:
            raise AttributeError("'h5f' must be a HDF5 file path or a valid h5py.File object. Got '%s'" % type(h5f))

        # Get base Group
        group = f[where]

        return f, group, close_when_done

    @classmethod
    def _h5_deserialize(cls, h5group, version):
        raise NotImplementedError()

    @classmethod
    def _h5_classname_read(cls, h5group):
        if "__class_name__" not in h5group.attrs:
            return None
        return h5group.attrs["__class_name__"]

    @classmethod
    def _h5_classname_write(cls, h5group):
        h5group.attrs["__class_name__"] = cls.__name__

    @classmethod
    def _h5_version_write(cls, h5group):
        h5group.attrs["__pymses_h5_version__"] = cls._pymses_h5_version

    @classmethod
    def _h5_version_read(cls, h5group):
        if "__pymses_h5_version__" not in h5group.attrs:
            return 1
        return h5group.attrs['__pymses_h5_version__']

    @classmethod
    def _fact_all_subclasses(cls, store_cache=True):
        """
        Returns a list of all sub(sub*)classes of a given class. Handles a subclass list cache dictionary.

        Parameters
        ----------
        store_cache; ``bool``
            Whether the list of subclasses of the class must be stored in a cache dictionary or not. Default True.

        Returns
        -------
        subclass_list: ``list``
            list of all subclasses of a given class
        """
        if cls in HDF5Serializable.__fact_subclasses_dict:
            return HDF5Serializable.__fact_subclasses_dict[cls]

        subclass_list = cls.__subclasses__() + [c for subclass in cls.__subclasses__()
                                                for c in subclass._fact_all_subclasses(store_cache=False)]

        if store_cache:
            HDF5Serializable.__fact_subclasses_dict[cls] = subclass_list
        return subclass_list

    @classmethod
    def from_HDF5(cls, h5fg, where="/"):
        """
        Returns an serialisable object from a HDF5 file.

        Parameters
        ----------
        h5fg: ``str`` or ``h5py.File`` or ``h5py.Group``
            HDF5 (h5py) file/group or filename.
        where: ``string``
            location of the object in the HDF5 tree. Default is tree root ("/").

        Raises
        ------
        imperr: ``ImportError``
            raises an ImportError if the 'h5py' module is not found.
        atterr: ``AttributeError``
            raises an AttributeError if any attribute is not valid.
        ioerr: ``IOError``
            raises an IOError if an error occured while reading the HDF5 file.
        """
        if not _h5py_available:
            raise ImportError("h5py module is required to handle HDF5 files in PyMSES !")

        if isinstance(h5fg, h5py.Group) and h5fg.__class__ == h5py.Group:
            group = h5fg
            close_when_done = False
            f = None
        else:
            f, group, close_when_done = cls.__open_h5file(h5fg, where)

        try:
            subclasses = cls._fact_all_subclasses()
            class_name = cls._h5_classname_read(group)
            # Instance class is a subclass of the current class => find the right class
            if class_name is not None and class_name != cls.__name__:
                o = None
                for kl in subclasses:
                    if kl.__name__ == class_name:
                        o = kl.from_HDF5(group)
                        break
                if o is None:
                    raise IOError("Cannot instantiate object of class '%s' in '%s'." % (class_name, group.name))
            else:
                # Read version number
                v = cls._h5_version_read(group)

                # Read object
                o = cls._h5_deserialize(group, v)
        except Exception as exc:
            raise IOError("HDF5 I/O error : %s" % exc)
        finally:
            if close_when_done and f is not None:
                f.close()

        return o


__all__ = ["HDF5Serializable"]
