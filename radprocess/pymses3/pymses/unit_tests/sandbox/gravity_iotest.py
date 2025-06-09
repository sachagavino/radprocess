#!/usr/bin/env python
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

import os
import  numpy  as N
import numpy
import argparse
from pymses.sources.ramses import _read_ramses


parser = argparse.ArgumentParser(description="Gravity I/O test script")
parser.add_argument("out_dir", help="RAMSES output repository")
parser.add_argument("--iout", type=int, help="Ramses output number (default 1)", dest="iout", metavar="iout", default=1)
parser.add_argument("--cpu", type=int, help="cpu file numer (default 1)", dest="icpu", metavar="icpu", default=1)
parser.add_argument("-s", help="Enable byteswapping", action="store_true", dest="swap")
args = parser.parse_args()

dir = os.path.abspath(args.out_dir)
odir = os.path.join(dir, "output_%5.5d" % args.iout)
grav_fname = os.path.join(odir, "grav_%5.5d.out%5.5d" % (args.iout, args.icpu))
hydro_fname = os.path.join(odir, "hydro_%5.5d.out%5.5d" % (args.iout, args.icpu))
amr_fname = os.path.join(odir, "amr_%5.5d.out%5.5d" % (args.iout, args.icpu))

print("Reading AMR file '%s'..." % amr_fname)
amrhdr, amrstruct = _read_ramses.read_amr(amr_fname, 99, int(args.swap))
print("amr_ngrids : %d" % amrstruct["ngrids"])
print("amr_readlmax : %d" % amrstruct["readlmax"])

print("Reading gravity file '%s'..." % grav_fname)
# Read [phi, gx, gy, gz] variables
data = _read_ramses.read_cells("grav", grav_fname, amrstruct["ndim"], amrstruct["ngrids"], amrstruct["readlmax"],
                               [0, 1, 2, 3], int(args.swap))

twotondim, ngrids = data.shape[1:]
print("ngrids = %d" % ngrids)
print("twotondim = %d" % twotondim)
dset = N.swapaxes(data, 0, 2).reshape(ngrids*twotondim, -1)
print(dset.shape)
# print dset[140:155, :]
print("Mean field : ", numpy.mean(dset, axis=0))

print("Reading hydro file '%s'..." % hydro_fname)
data = _read_ramses.read_cells("hydro", hydro_fname, amrstruct["ndim"], amrstruct["ngrids"], amrstruct["readlmax"],
                               [0, 1, 2, 3], int(args.swap))

twotondim, ngrids = data.shape[1:]
print("ngrids = %d" % ngrids)
print("twotondim = %d" % twotondim)
dset = N.swapaxes(data, 0, 2).reshape(ngrids*twotondim, -1)
print(dset.shape)
# print dset[140:155, :]
mask = dset[:, 0]>0.0
print("Mean field : ", numpy.mean(dset[mask, 0], axis=0))
print("min/max : ", numpy.min(dset[mask, 0]), numpy.max(dset[mask, 0]))
print("std. deviation : ", numpy.std(dset[mask, 0], axis=0))

