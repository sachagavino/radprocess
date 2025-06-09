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

from numpy import array
import os
from pymses.analysis.visualization import raytracing, MaxLevelOperator, Camera, image_plot_utils as ImgPlot
from pymses import RamsesOutput
from pymses.utils import constants as C

# Ramses data
ioutput = 193
datadir = os.path.join(os.path.expanduser("~"), "data/Aquarius/output")
ro = RamsesOutput(datadir, ioutput)

# Map operator : max. AMR level of refinement along the line-of-sight
scal_op = MaxLevelOperator()

# Map region
center = [0.567811, 0.586055, 0.559156]
axes = {"los": array([-0.172935, 0.977948, -0.117099])}

# Map processing
rt = raytracing.RayTracer(ro, ["rho"])
for axname, axis in axes.items():
    cam = Camera(center=center, line_of_sight_axis=axis, up_vector="z", region_size=[4.0E-2, 4.0E-2],
                 size_unit=ro.info["unit_length"], distance=2.0E-2, far_cut_depth=2.0E-2, map_max_size=512,
                 log_sensitive=False)
    map = rt.process(scal_op, cam)
    scale = ro.info["unit_length"].express(C.Mpc)

    # Save map into HDF5 file
    mapname = "gas_rt_lmax_%s_%5.5i" % (axname, ioutput)
    h5fname = ImgPlot.save_map_HDF5(map, cam, map_name=mapname)

    # Plot map into Matplotlib figure/PIL Image
    fig = ImgPlot.save_HDF5_to_plot(h5fname, map_unit=("AMR level", 1.0), axis_unit=("Mpc", scale), cmap="jet",
                                    discrete=True)
    # pil_img = save_HDF5_to_img(h5fname, cmap="jet", discrete=True)

# Save into PNG image file
# save_HDF5_to_plot(h5fname, map_unit=("AMR level",1.0), axis_unit=("Mpc", scale), img_path="./", cmap="jet", discrete=True)
# save_HDF5_to_img(h5fname, img_path="./", cmap="jet", discrete=True)

# pylab.show()
