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
:mod:`pymses.sources.ramses.sink` --- Ramses sink particles module
------------------------------------------------------------------

"""
from . import filename_utils
import numpy as N
from pymses.core.datasets import PointDataset
from pymses.etc.rcsetup import rcConfiguration


class SinkParticleDataReader(object):
    """
    Sink particle data reader class

    Parameters
    ----------
    output_path: ``string``
        Ramses output base directory path.
    ioutput: ``int``
        Ramses output number.
    """
    def __init__(self, output_path, ioutput):
        super(SinkParticleDataReader, self).__init__()

        self._sink_filename = filename_utils.sink_filename(output_path, ioutput, check_exists=True)
        self._sink_data = None

    @property
    def sink_filepath(self):
        return self._sink_filename

    def read(self, sink_fields):
        """
        TODO
        :param sink_fields:
        :return:
        """
        # Make sure the sink particle position field is read
        cfg = rcConfiguration().Ramses.sink_fields
        read_fields = sink_fields[:]
        if cfg.sink_position_field not in read_fields:
            read_fields.append(cfg.sink_position_field)

        ivars, flist = cfg.gather_read_fields(read_fields)

        # Read sink particle file and init PointDataset
        d = N.loadtxt(self._sink_filename, skiprows=5, ndmin=2, comments=" ==", usecols=set(ivars))
        pdset = None
        for field in flist:
            if field.field_name == cfg.sink_position_field:
                pdset = PointDataset(field.collect_data(d))
                break

        # Security check
        if pdset is None:  # Should never happen
            raise IOError("Cannot read sink particle position field from"
                          " '{filepath!s}'.".format(filepath=self._sink_filename))

        # Add fields to the PointDataset
        for field in flist:
            if field.field_name != cfg.sink_position_field:
                field.gather(d, pdset)

        return pdset


__all__ = ['SinkParticleDataReader']
