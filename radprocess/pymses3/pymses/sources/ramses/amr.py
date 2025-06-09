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
:mod:`pymses.sources.ramses.amr` --- flexible raw RAMSES AMR file reading module
--------------------------------------------------------------------------------
"""
from . import _read_ramses


def read_ramses_amr_file(amr_filename, max_read_level=None, swap=False):
    """ Reads a RAMSES .amr file into memory

    Arguments:

        amr_filename -- filename of the .amr file

        max_read_level -- read all levels <= max_read_level (default: read all)

        swap -- swap bytes (to match different endianness)

    Returns: (header_dict, amr_struct_dict)
    """
    if max_read_level is None:
        max_read_level = 99

    try:
        out = _read_ramses.read_amr(amr_filename, max_read_level, int(swap))
    except _read_ramses.RamsesIOError as err:
        raise IOError(str(err))
    return out
