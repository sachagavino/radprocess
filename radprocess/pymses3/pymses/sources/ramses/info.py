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
info.py -- RAMSES .info reading routine
"""

import re
import os
import ast
import numpy
import pymses.utils.constants as C
from ._read_ramses import read_hydro_header


def read_ramses_info_file(info_filename):
    """ Reads a RAMSES .info file, returning a dictionary of parameter values
    """
    if not os.path.isfile(info_filename):
        raise AttributeError("Cannot find RAMSES output info file '%s'." % info_filename)

    re_kwline = re.compile(r"^\w+\s*=\s*.+$")
    re_split = re.compile(r"\s*=\s*")
    re_split2 = re.compile(r"\s+")
    re_emptyline = re.compile(r"^\s*$")
    re_header_boundkeys = re.compile(r"^\s*DOMAIN\s*ind_min\s*ind_max\s*$")

    bound_table = None
    par_dict = {}
    is_param = True

    with open(info_filename, 'r') as info_fileobj:
        # Info.txt file parsing
        for line in info_fileobj.readlines():
            if re_emptyline.match(line) is not None:
                continue
            if re_header_boundkeys.match(line) is not None:
                is_param = False
                continue
            if re_kwline.match(line) is not None and not is_param:
                is_param = True
            if is_param:
                # Parameter parsing
                params = re_split.split(line)
                par_name = params[0].strip()
                par_val = params[1].strip().strip("\n")
                if par_name == "ordering type":
                    par_dict["ordering"] = par_val
                    bound_table = []
                else:
                    par_dict[par_name] = par_val
            else:
                # Boundary keys parsing
                params = re_split2.split(line)
                bound_table.append([float(params[2]), float(params[3].strip().strip("\n"))])

    info_dict = {}
    for key, value in par_dict.items():
        if key in ["ncpu", "ndim", "levelmin", "levelmax", "ngridmax", "nstep_coarse"]:
            info_dict[key] = int(value)
        elif key in ["boxlen", "time", "aexp", "H0", "omega_m", "omega_l", "omega_k", "omega_b"]:
            info_dict[key] = float(value)
        elif key == "ordering type":
            info_dict["ordering"] = value
        elif key == "unit_l":
            unit_l = float(value)
            info_dict["unit_length"] = unit_l * C.cm
        elif key == "unit_d":
            unit_d = float(value)
            info_dict["unit_density"] = unit_d * C.g / C.cm ** 3
        elif key == "unit_t":
            unit_t = float(value)
            info_dict["unit_time"] = unit_t * C.s
        else:  # Custom (user's RAMSES patch) parameters
            try:
                info_dict[key] = ast.literal_eval(value)
            except ValueError:  # Cannot evaluate value type, maybe it really was a string ?
                info_dict[key] = value

    # -------------------------------------- Unit post-processing --------------------------------------------- #
    info_dict["unit_velocity"] = info_dict["unit_length"] / info_dict["unit_time"]
    info_dict["unit_pressure"] = info_dict["unit_density"] * info_dict["unit_velocity"] ** 2
    info_dict["unit_temperature"] = info_dict["unit_velocity"] ** 2 * C.mH / C.kB
    info_dict["unit_mag"] = numpy.sqrt(4 * numpy.pi * unit_d * (unit_l / unit_t) ** 2) * C.Gauss
    info_dict["unit_gravpot"] = info_dict["unit_velocity"]**2
    # To use only for collisionless particles (needs to ben done BEFORE 'unit_length' renormalisation by boxlen).
    # Bug introduced by revision 897 (May, 16, 2013)
    # see issues #6 (https://bitbucket.org/dchapon/pymses/issues/6/)
    #            #10 (https://bitbucket.org/dchapon/pymses/issues/10/).
    info_dict["unit_mass"] = info_dict["unit_density"] * info_dict["unit_length"] ** 3

    # Boxlen renormalisation
    info_dict["unit_length"] = info_dict["unit_length"] * info_dict["boxlen"]
    info_dict["unit_gravpot"] = info_dict["unit_gravpot"] * info_dict["boxlen"]
    # --------------------------------------------------------------------------------------------------------- #

    if info_dict["ordering"] == "hilbert" or info_dict["ordering"] == "planar":
        keys = numpy.zeros(info_dict["ncpu"] + 1)
        bound_table = numpy.asarray(bound_table)
        keys[0] = bound_table[0, 0]
        keys[1:] = bound_table[:, 1]
        if info_dict["ordering"] == "hilbert":
            info_dict["dom_decomp_Hilbert_keys"] = keys
        else:  # planar
            info_dict["dom_decomp_planar_keys"] = keys
            info_dict["dom_decomp"] = None
    else:
        info_dict["dom_decomp"] = None

    return info_dict


def read_ramses_hydro_header(hydro_filename, swap=False):
    try:
        hydro_hdr = read_hydro_header(hydro_filename, int(swap))
    except _read_ramses.RamsesIOError as err:
        raise IOError(str(err))

    return {"nvar_file": hydro_hdr["nvar_file"],
            "gamma": hydro_hdr["gamma"]}


__all__ = ["read_ramses_info_file", "read_ramses_hydro_header"]
