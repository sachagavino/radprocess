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
:mod:`pymses.utils.filename`: file path utility functions module
----------------------------------------------------------------
"""
import os
import re


class FileUtil(object):
    """
    File path utility abstract class
    """
    HDF5_FILE = 1
    PNG_FILE = 2
    FITS_FILE = 3

    @classmethod
    def valid_filepath(cls, filename, append_extension=None):

        if not isinstance(filename, str):
            raise filename("'fname' attribute must be a file path (string). Got %s." % type(filename))

        # Get file absolute path
        fname = os.path.abspath(filename)

        if append_extension is None:
            return fname
        elif append_extension == FileUtil.HDF5_FILE:
            h5_file_regex = re.compile("^.+\.(h5|H5||hdf5|HDF5)$")
            if h5_file_regex.match(fname) is None:  # Append .h5 extension to filename
                fname = "%s.h5" % fname
        elif append_extension == FileUtil.PNG_FILE:
            png_file_regex = re.compile("^.+\.(png|PNG)$")
            if png_file_regex.match(fname) is None:  # Append .png extension to filename
                fname = "%s.png" % fname
        elif append_extension == FileUtil.FITS_FILE:
            fits_file_regex = re.compile("^.+\.(fits|FITS)$")
            if fits_file_regex.match(fname) is None:  # Append .fits extension to filename
                fname = "%s.fits" % fname
        else:
            print("Warning: unhandled file extension ('%s')" % append_extension)
        return fname

    @classmethod
    def new_filepath(cls, fname, append_extension=None, create_dir=True):
        valid_fname = cls.valid_filepath(fname, append_extension)

        file_dir = os.path.dirname(valid_fname)
        if create_dir and not os.path.isdir(file_dir):  # Create directory if it does not exist
            os.mkdir(file_dir)
        return valid_fname


__all__ = ["FileUtil"]
