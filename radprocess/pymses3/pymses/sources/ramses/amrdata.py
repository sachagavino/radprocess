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
:mod:`pymses.sources.ramses.amrdata` --- flexible raw RAMSES AMR data file reading module
-----------------------------------------------------------------------------------------

"""
from ._read_ramses import read_cells, RamsesIOError
from numpy import swapaxes, ascontiguousarray


def read_ramses_amr_data_file(data_filename, is_hydro, amr, ivars_to_read, swap=False, grav_compat=False):
    r"""
    Reads a RAMSES AMR data file into memory (.hydro, .grav, ...)

    Parameters
    ----------
    data_filename : ``string``
        filename of the data file.
    is_hydro      : ``int``
        0 if data file type is "grav", otherwise 1.
    amr           : AMR structure
        the AMR structure, as output by :func:`~pymses.sources.ranses.amr.read_ramses_amr_file`.
    ivars_to_read : ``list`` of ``int``
        list of variable ids to read.
    swap          : ``bool``
        swap bytes if True (to match different endianness).
    grav_compat: ``bool``
        Old Ramses versions (prior to commit 5dd90f3, 2012-10-04) and new Ramses versions (posterior to commit bce4454,
        2015-07-09) work with nvar_file (number of scalar fields written in the file) parameter written in file header :
            \_,-> gravitational attraction vector field only (old version)
            \_,-> potential scalar field (phi) + gravitational attraction vector field (new version)

        For Ramses versions BETWEEN 2012-10-04 AND 2015-07-09, ndim + 1 variables were saved but only ndim was written
        in the file header.
        To read output files written by this flavor of Ramses, set the `grav_compat` attribute to True. (Default False).

    Returns
    -------
    a : ``numpy.ndarray``
        a C-ordered, contiguous Numpy ndarray of shape (ngrids, twotondim, nvars) containing the amr data values for the
        required fields.

    """
    dummy_amrhdr, amrstruct = amr

    try:
        data = read_cells(data_filename, is_hydro, amrstruct["ndim"], amrstruct["ngrids"],
                          amrstruct["readlmax"], ivars_to_read, int(swap), int(grav_compat))
    except RamsesIOError as exc:
        raise IOError(exc)

    # Low-level C I/O routines assemble data as a contiguous, C-ordered (nvars, twotondim, ngrids) numpy.ndarray
    # Swap data => shape : (ngrids, twotondim, nvars)
    ####### WARNING : must keep C-ordered contiguity !!! #######
    return ascontiguousarray(swapaxes(data, 0, 2))
