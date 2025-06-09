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
:mod:`pymses.etc.rcsetup` --- Runtime changes configuration setup module
------------------------------------------------------------------------

the ~/.pymses/pymsesrc file can be edited to override the default PyMSES configuration parameters. This file is
written in the JSON format

Example
-------

To set seven custom chemical concentrations 'Xi' fields in any AMR datasource read by PyMSES, you may edit your
~/.pymses/pymsesrc file as follows :
{
    "Version": 1,
    "Multiprocessing max. nproc": 8,
    "RAMSES":{
        "ndimensions": 3,
        "amr_field_descr": [
            {"__type__": "scalar_field", "__file_type__": "hydro", "name": "rho", "ivar": 0},
            {"__type__": "vector_field", "__file_type__": "hydro", "name": "vel", "ivars": [1, 2, 3]},
            {"__type__": "vector_field", "__file_type__": "hydro", "name": "Bl", "ivars": [4, 5, 6]},
            {"__type__": "vector_field", "__file_type__": "hydro", "name": "Br", "ivars": [7, 8, 9]},
            {"__type__": "scalar_field", "__file_type__": "hydro", "name": "P", "ivar": 10},
            {"__type__": "multivalued_field", "__file_type__": "hydro", "name": "Xi", "ivar_first": 11, "nb_vars": 7},
            {"__type__": "scalar_field", "__file_type__": "grav", "name": "phi", "ivar": 0},
            {"__type__": "vector_field", "__file_type__": "grav", "name": "g", "ivars": [1, 2, 3]}
        ]
    }
}

"""
from pymses.sources.field_types import Scalar, Vector, MultiValued
from pymses.sources.ramses.field_descr import AmrDataFileFieldDescriptor, SinkParticleFileFieldDescriptor,\
    RamsesAmrScalar, RamsesAmrVector, RamsesAmrMultiValued
import os
import numpy as N
# from shutil import copyfile
import json


def _validate_rcconfig(base_cfg, custom_cfg):
    """
    Update a base configuration dictionary with values from a custom configuration dictionary
    """
    keys = list(base_cfg.keys())
    for k in custom_cfg:
        if k not in keys:  # Add custom (key, value) pair in default dict
            base_cfg[k] = custom_cfg[k]
        else:  # Same key found in both dict objects
            if type(base_cfg[k]) is dict:
                if type(custom_cfg[k]) is dict:  # recursive call
                    _validate_rcconfig(base_cfg[k], custom_cfg[k])
                else:  # Ignore user value if it is not a dict
                    pass
            elif isinstance(base_cfg[k], RamsesConfiguration):  # Updates Ramses configuration
                base_cfg[k].update(custom_cfg[k])
            else:  # Overriding value in default dict for key k
                base_cfg[k] = custom_cfg[k]


class rcConfiguration(object):
    __instance = None
    _version_key = "Version"
    _multiprocessing_max_nproc_key = "Multiprocessing max. nproc"
    _ramses_key = "RAMSES"
    _pymsesrc_fname = "pymsesrc"
    _dot_pymses_dir = ".pymses"

    def __new__(cls):
        """
        Load the PyMSES runtime changes configuration file. Get the location of the file (JSON format) as follows :
        * `$PWD/pymsesrc`, if it exists
        * `$HOME/.pymses/pymserc`, if it exists
        * Lastly, it takes the `pymses/etc/pymsesrc` as the default configuration file.
        """
        if rcConfiguration.__instance is not None:
            return rcConfiguration.__instance

        rcConfiguration.__instance = super(rcConfiguration, cls).__new__(cls)

        # Default pymsesrc configuration file
        _etc_pymsesrc = rcConfiguration._get_etc_pymsesrc_filepath()

        # Local directory's `pymsesrc` file
        _fname = rcConfiguration._get_local_pymsesrc_filepath()
        if not os.path.isfile(_fname):
            # Look for a '~/.pymses/pymsesrc' file in the user's home directory
            _fname = rcConfiguration._get_user_dot_pymses_pymsesrc_filepath()
            if not os.path.isfile(_fname):
                # Using default `pymses/etc/pymsesrc` configuration file
                # copyfile(_etc_pymsesrc, _fname)
                _fname = None

        rcConfiguration.__instance.__load_json_cfg_file(_fname)
        return rcConfiguration.__instance

    def __load_json_cfg_file(self, filename=None):

        self._user_pymsesrc_dict = {}
        self._pymsesrc_dict = {}
        self._ramses_config = None

        def as_pymsesrc_dict(d):
            if rcConfiguration._multiprocessing_max_nproc_key in d:
                if not isinstance(d[rcConfiguration._multiprocessing_max_nproc_key], int):
                    raise TypeError("'rcConfiguration._multiprocessing_max_nproc_key' value must be an int")
            if rcConfiguration._ramses_key in d:
                self._ramses_config = RamsesConfiguration(d[rcConfiguration._ramses_key])
                del d[rcConfiguration._ramses_key]
                return d
            # elif rcConfiguration._other_top_level_config_key in d:
            #     pass
            else:
                return d

        with open(rcConfiguration._get_etc_pymsesrc_filepath()) as pymsesrc_dict_file:
            self._pymsesrc_dict = json.load(pymsesrc_dict_file, object_hook=as_pymsesrc_dict)

        if filename is not None:
            with open(filename) as default_pymsesrc_dict_file:
                try:
                    self._user_pymsesrc_dict = json.load(default_pymsesrc_dict_file, object_hook=as_pymsesrc_dict)
                except:
                    print("Warning : invalid JSON file format. '%s' file will be ignored and default rc settings " \
                          "will be used. " % filename)

            _validate_rcconfig(self._pymsesrc_dict, self._user_pymsesrc_dict)

    @classmethod
    def _get_etc_pymsesrc_filepath(cls):
        """
        Get the system default pymsesrc config file path

        :return: ``string``
            system pymsesrc file path
        """
        _etc_pymsesrc = os.path.join(os.path.dirname(__file__), rcConfiguration._pymsesrc_fname)
        return _etc_pymsesrc

    @classmethod
    def _get_user_dot_pymses_pymsesrc_filepath(cls):
        """
        Get the user ~/.pymses/pymsesrc config file path

        :return: ``string``
            user pymsesrc file path
        """
        # Look for a '.pymses/' directory in the user's home directory
        _home_dir = os.path.expanduser("~")
        _dot_pymses_path = os.path.join(_home_dir, rcConfiguration._dot_pymses_dir)
        if not os.path.isdir(_dot_pymses_path):
            print("Creating '%s' directory into the home directory : '%s'." % (_dot_pymses_path, _home_dir))
            os.makedirs(_dot_pymses_path)  # Create `~/.pymses` directory
        _user_pymsesrc = os.path.join(_dot_pymses_path, rcConfiguration._pymsesrc_fname)
        return _user_pymsesrc

    @classmethod
    def _get_local_pymsesrc_filepath(cls):
        """
        Get the local $PWD/pymsesrc config file path

        :return: ``string``
            local pymsesrc file path
        """
        #  Local directory's `pymsesrc` file
        _local_pymsesrc = os.path.join(os.getcwd(), rcConfiguration._pymsesrc_fname)
        return _local_pymsesrc

    @property
    def Ramses(self):
        """
        Get the Ramses run-level changes configuration settings (:class:`~pymses.etc.rcsetup.RamsesConfiguration``
        instance)
        """
        return self._ramses_config

    @property
    def multiprocessing_max_nproc(self):
        return self._pymsesrc_dict[rcConfiguration._multiprocessing_max_nproc_key]

    @multiprocessing_max_nproc.setter
    def multiprocessing_max_nproc(self, nproc_max):
        if not isinstance(nproc_max, int) or nproc_max < 1:
            raise ValueError("Invalid 'multiprocessing_max_nproc' value.")
        self._pymsesrc_dict[rcConfiguration._multiprocessing_max_nproc_key] = nproc_max


class RamsesConfiguration(object):
    _amr_field_descr_key = "amr_field_descr"
    _sink_field_descr = "sink_field_descr"
    _ndimensions_key = "ndimensions"
    _field_type_key = "__type__"
    _datatype = "__dtype__"
    _sink_position_flag = "position"
    _filetype_key = "__file_type__"

    def __init__(self, ramses_dict):
        self.ramses_config_dict = ramses_dict

        # AMR field descriptor
        amr_field_list = self.ramses_config_dict.pop(RamsesConfiguration._amr_field_descr_key, None)
        self._amr_field_descr = None
        self._init_amr_field_descriptor(amr_field_list)

        # Sink particles field descriptor
        sink_field_list = self.ramses_config_dict.pop(RamsesConfiguration._sink_field_descr, None)
        self._sink_field_descr = None
        self._init_sink_field_descriptor(sink_field_list)

    def _init_amr_field_descriptor(self, amr_field_list=None):
        if amr_field_list is None:
            return

        # Build RamsesAmr* field list
        flist = []
        for field_dict in amr_field_list:
            if 'name' not in field_dict:
                raise IOError("Missing AMR data field 'name' value.")

            fname = field_dict["name"]

            if RamsesConfiguration._field_type_key not in field_dict:
                errMsg = "Missing field type '{field_type_key!s}' value for AMR data field " \
                         "named '{field_name!s}'.".format(field_type_key=RamsesConfiguration._field_type_key,
                                                          field_name=fname)
                raise IOError(errMsg)

            t = field_dict[RamsesConfiguration._field_type_key]

            if RamsesConfiguration._filetype_key not in field_dict:
                errMsg = "Missing file type '{file_type_key!s}' value for AMR data field " \
                         "named '{field_name!s}'.".format(file_type_key=RamsesConfiguration._filetype_key,
                                                          field_name=fname)
                raise IOError(errMsg)

            pfiletype = field_dict[RamsesConfiguration._filetype_key]

            if t == RamsesAmrScalar.field_type():  # Scalar AMR field
                if 'ivar' not in field_dict:
                    raise IOError("Missing AMR data scalar field 'ivar' integer value "
                                  "for field named {field_name!s}.".format(field_name=fname))

                flist.append(RamsesAmrScalar(pfiletype, fname, field_dict['ivar']))
            elif t == RamsesAmrVector.field_type():  # Vector AMR field
                if 'ivars' not in field_dict:
                    raise IOError("Missing AMR data scalar field 'ivars' integer value list "
                                  "for field named {field_name!s}.".format(field_name=fname))

                flist.append(RamsesAmrVector(pfiletype, fname, field_dict['ivars']))
            elif t == RamsesAmrMultiValued.field_type():  # Multivalued AMR field
                if 'ivar_first' not in field_dict:
                    raise IOError("Missing AMR data multivalued field 'ivar_first' integer value "
                                  "for field named {field_name!s}.".format(field_name=fname))
                if 'nb_vars' not in field_dict:
                    raise IOError("Missing AMR data multivalued field 'nb_vars' integer value "
                                  "for field named {field_name!s}.".format(field_name=fname))

                flist.append(RamsesAmrMultiValued(pfiletype, fname, field_dict['ivar_first'], field_dict['nb_vars']))

        self._amr_field_descr = AmrDataFileFieldDescriptor(flist)

    def _init_sink_field_descriptor(self, sink_field_list=None):
        if sink_field_list is None:
            return

        # Build sink particle Scalar/Vector/MultiValued field list
        sflist = []
        pos_field = None
        for field_dict in sink_field_list:
            if 'name' not in field_dict:
                raise IOError("Missing sink particle data field 'name' value.")

            fname = field_dict["name"]

            if RamsesConfiguration._field_type_key not in field_dict:
                errMsg = "Missing field type '{field_type_key!s}' value for sink particle data field " \
                         "named '{field_name!s}'.".format(field_type_key=RamsesConfiguration._field_type_key,
                                                          field_name=fname)
                raise IOError(errMsg)

            t = field_dict[RamsesConfiguration._field_type_key]

            if RamsesConfiguration._datatype not in field_dict:
                errMsg = "Missing datatype '{dtype_key!s}' value for sink particle data field " \
                         "named '{field_name!s}'.".format(dtype_key=RamsesConfiguration._datatype,
                                                          field_name=fname)
                raise IOError(errMsg)

            dt = field_dict[RamsesConfiguration._datatype]
            if dt == 'int':
                datatype = N.int32
            elif dt == 'float':
                datatype = N.float32
            elif dt == 'double':
                datatype = N.float64
            elif dt == 'bool':
                datatype = bool
            else:
                raise IOError("Invalid datatype '{dtype!s}' for sink particle data field named '{field_name!s}'."
                              "Valid datatypes are 'int', 'float', 'double' and 'bool'.".format(dtype=dt,
                                                                                                field_name=fname))

            # Get sink particle position field name
            is_pos_field = field_dict.get(RamsesConfiguration._sink_position_flag, False)
            if is_pos_field:
                if pos_field is not None:  # Assert sink particle position field unicity
                    raise IOError("Multiple sink particle position fields definition " \
                                  "'{pos_field1!s}' and '{pos_field2!s}'".format(pos_field1=pos_field,
                                                                                 pos_field2=fname))

                pos_field = fname

            if t == Scalar.field_type():  # Scalar field
                if 'ivar' not in field_dict:
                    raise IOError("Missing sink particle data scalar field 'ivar' integer value "
                                  "for field named {field_name!s}.".format(field_name=fname))

                sflist.append(Scalar(fname, field_dict['ivar'], datatype))
            elif t == Vector.field_type():  # Vector field
                if 'ivars' not in field_dict:
                    raise IOError("Missing sink particle data vector field 'ivars' integer value list "
                                  "for field named {field_name!s}.".format(field_name=fname))

                sflist.append(Vector(fname, field_dict['ivars'], datatype))
            elif t == MultiValued.field_type():  # Multivalued field
                if 'ivar_first' not in field_dict:
                    raise IOError("Missing sink particle data multivalued field 'ivar_first' integer value "
                                  "for field named {field_name!s}.".format(field_name=fname))
                if 'nb_vars' not in field_dict:
                    raise IOError("Missing sink particle data multivalued field 'nb_vars' integer value "
                                  "for field named {field_name!s}.".format(field_name=fname))

                sflist.append(MultiValued(fname, field_dict['ivar_first'], field_dict['nb_vars'], datatype))

        if pos_field is None:
            raise IOError("No sink particle position field defined.")

        self._sink_field_descr = SinkParticleFileFieldDescriptor(sflist, pos_field)

    @property
    def ndims(self):
        """
        Get the number of space dimensions of the Ramses simulation output data

        :rtype : int
        :return: the number of space dimensions within the Ramses output data (or None, if the parameter is not defined)
        """
        return self.ramses_config_dict.get(RamsesConfiguration._ndimensions_key)

    @property
    def amr_fields(self):
        """
        Get the AMR data field descriptor

        :rtype : :class:`~pymses.sources.ramses.field_descr.AmrDataFileFieldDescriptor`
        :return: the AMR data field list  descriptor
        """
        return self._amr_field_descr

    @property
    def sink_fields(self):
        """
        Get the sink particle data field descriptor

        :rtype : :class:`~pymses.sources.ramses.field_descr.SinkParticleFileFieldDescriptor`
        :return: the sink particle data field list descriptor
        """
        return self._sink_field_descr

    def update(self, user_cfg):
        """
        Updates the current :class:`~pymses.etc.rcsetup.RamsesConfiguration`` object with the parameters of a given
        user-defined :class:`~pymses.etc.rcsetup.RamsesConfiguration`` instance

        :param user_cfg: user-defined :class:`~pymses.etc.rcsetup.RamsesConfiguration`` instance
        """
        _validate_rcconfig(self.ramses_config_dict, user_cfg.ramses_config_dict)
        if not isinstance(user_cfg, RamsesConfiguration):
            return

        # Overrides the AMR field list + sink particle data field list
        self._amr_field_descr = user_cfg.amr_fields
        self._sink_field_descr = user_cfg.sink_fields


__all__ = ["rcConfiguration"]
