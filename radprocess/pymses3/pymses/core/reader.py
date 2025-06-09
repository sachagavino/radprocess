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
:mod:`pymses.core.reader` --- PyMSES abstract multiple data file reader module
==============================================================================

"""


class MultiFileDataReader(object):
    """
    Abstract multiple-file data reader class
    """
    def read_file(self, idata, max_res):
        """
        Abstract reader data file reading method

        Parameters
        ----------
        idata : ``int``
            data file index
        max_res : ``int``
            maximum spatial resolution to read from the data file

        """
        raise NotImplementedError()


__all__ = ["MultiFileDataReader"]
