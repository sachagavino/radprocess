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

from PIL import Image
import os

if __name__ == "__main__":
    dirname = "./_build/plot_directive/pyplots"
    for file in os.listdir(dirname):
        (file_name, file_ext) = os.path.splitext(file)
        if file_ext == '.png':
            (file_name, file_ext) = os.path.splitext(file_name)
            if file_ext == '.hires':
                im = Image.open("%s/%s" % (dirname,file))
                rgb_im = im.convert("RGB")
                rgb_im.save("%s/%s.pdf" % (dirname,file_name))
                print("%s/%s.pdf created" % (dirname,file_name))
