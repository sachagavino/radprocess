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

# Colormap map testing script
from argparse import ArgumentParser

from pymses.analysis import DataMap
from pymses.analysis.plot import Plot2D, PlainPNGImage

parser = ArgumentParser(description="Apply a colormap to a Datamap HDF5 file and save it into an image or a plot")
parser.add_argument('h5_fname', help="Datamap HDF5 file name")
parser.add_argument('out_fname', help="Output file name")
parser.add_argument('-p', '--plot', action='store_true', help='Save a plot instead of an image')
parser.add_argument('-c', dest='cmap_name', default='BlackBlueYellowRed', metavar='cmap_name',
                    help="required colormap name (default 'BlackBlueYellowRed'")

# Parse arguments
args = parser.parse_args()
HDF5fileName = args.h5_fname
output_filename = args.out_fname
cmap_name = args.cmap_name

# Load DataMap
dmap = DataMap.from_HDF5(HDF5fileName)

# Plot
if args.plot:
    plot_gen = Plot2D(cmap=cmap_name)
    dmap.save_plot(plot_gen, plot_fname=output_filename, fraction=1.0)
else:
    img_gen = PlainPNGImage(cmap=cmap_name)
    dmap.save_PNG(img_gen, img_fname=output_filename, fraction=1.0)
