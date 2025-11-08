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
:mod:`pymses.sources.ramses.particles` --- utility functions for RAMSES outputs filename
----------------------------------------------------------------------------------------
"""
import os
import re

HYDRO_FILE_PREFIX = "hydro"
GRAV_FILE_PREFIX = "grav"
PART_FILE_PREFIX = "part"
AMR_FILE_PREFIX = "amr"
SINK_FILE_PREFIX = "sink"

SUPPORTED_FIELD_FILE_LIST = [HYDRO_FILE_PREFIX, GRAV_FILE_PREFIX, PART_FILE_PREFIX]


def output_path(output_repos, iout, check_exists=True):
    """
    Get a RAMSES output path from an output repository and output number

    Paramters
    ---------
    output_repos:
    iout:
    check_exists: ``bool``
        output directory existence needs to be checked ?

    Returns
    -------
    path; ``string``
        the complete path of the Ramses output directory
    """
    path = os.path.join(output_repos, "output_{ioutput:05d}".format(ioutput=iout))
    if check_exists and not os.path.isdir(path):
        raise IOError("Ramses output path '{routput_path!s}' does not exist.".format(routput_path=path))
    return path


def info_filename(output_repos, iout, check_exists=True):
    """
    Get a RAMSES output info_XXXXX.txt file name

    Parameters
    ----------
    output_repos:
    iout:
    check_exists: ``bool``
        info file existence needs to be checked ?

    Returns
    -------
    fname: ``string``
        the complete path of the Ramses output info file
    """

    path = output_path(output_repos, iout, check_exists)
    fname = os.path.join(path, "info_{ioutput:05d}.txt".format(ioutput=iout))
    if check_exists and not os.path.isfile(fname):
        raise IOError("Ramses output info file '{routput_info_fname!s}' does not exist.".format(routput_info_fname=fname))
    return fname


def namelist_filename_list(output_repos):
    """
    Search for namelist (*.nml) files within the output repository.

    Parameters
    ----------
    output_repos: output directory path

    Returns
    -------
    l: ``list`` or None.
        list of the namelist filenames found in the output repository or None if none were found.
    """
    if not os.path.isdir(output_repos):
        raise IOError("'{path!s}' is not a valid directory path.".format(path=output_repos))

    ls = os.listdir(output_repos)
    namelist_re = re.compile("^.+\.(nml|NML)$")
    l = []
    for file in ls:
        if not os.path.isfile(os.path.join(output_repos, file)):
            continue
        namelist_re.match(file)
        l.append(file)

    if len(l) == 0:
        return None
    else:
        return l


def sink_filename(output_repos, iout, check_exists=True):
    """
    Get a Ramses output sink file name

    Parameters
    ----------
    output_repos : ``string``
        the directory path containing the "output_*" directories
    iout : ``int``
        the output number
    check_exists: ``bool``
        sink file existence needs to be checked ?

    Returns
    -------
    sink_fname : ``string``
        the complete path of the Ramses output sink file
    """

    path = output_path(output_repos, iout, check_exists)
    sink_fname = os.path.join(path, "sink_{ioutput:05d}.out".format(ioutput=iout))
    if check_exists and not os.path.isfile(sink_fname):
        raise IOError("Ramses output sink particle file '{routput_sink_fname!s}' does not exist.".format(routput_sink_fname=sink_fname))
    return sink_fname


def amrlike_filename(ftype, output_repos, iout, icpu, check_exists=True):
    """
    Get a RAMSES per-CPU amr data file name.

    Parameters
    ----------
    ftype : ``string``
        the file type: "amr", "hydro", "grav", "part", etc
    output_repos : ``string``
        the directory path containing the "output_*" directories
    iout : ``int``
        the output number
    icpu : ``int``
        the CPU number
    check_exists: ``bool``
        amr data file existence needs to be checked ?

    Returns
    -------
    filename : ``string``
        the complete path of the amr data cpu file
    """
    path = output_path(output_repos, iout, check_exists)
    filename = os.path.join(path, "{fname_prefix!s}_{ioutput:05d}.out{cpu_number:05d}".format(fname_prefix=ftype,
                                                                                              ioutput=iout,
                                                                                              cpu_number=icpu))
    if check_exists and not os.path.isfile(filename):
        raise IOError("Ramses output amr data file '{amr_filename!s}' does not exist.".format(amr_filename=filename))
    return filename


def search_valid_outputs(out_dir):
    """
    Computes the ``int`` ``list`` of output number available in a given `out_dir`
    RAMSES outputs directory

    Parameters
    ----------
    out_dir : ``string``
        path of the directory containing all RAMSES outputs

    Returns
    -------
    ilist : ``list``
        sorted number list of the available outputs

    """
    ilist = []
    outdir_regexp = re.compile("output_[0-9]{5}")
    iout_regexp = re.compile("[0-9]{5}")
    ls = os.listdir(out_dir)
    for file in ls:
        res = outdir_regexp.findall(file)
        if len(res) > 0:
            if res[0] == file:
                ilist.append(int(iout_regexp.findall(file)[0]))

    ilist.sort()
    return ilist


__all__ = ["search_valid_outputs", "output_path", "info_filename", "amrlike_filename", "sink_filename",
           'namelist_filename_list', "HYDRO_FILE_PREFIX", "PART_FILE_PREFIX", "AMR_FILE_PREFIX", "SINK_FILE_PREFIX"]
