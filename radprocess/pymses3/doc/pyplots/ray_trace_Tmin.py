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

from numpy import array, zeros_like
import pylab
import os
from pymses.analysis.visualization import Operator, Camera, raytracing, image_plot_utils as ImgPlot
from pymses import RamsesOutput
from pymses.utils import constants as C

# Ramses data
ioutput = 193
datadir = os.path.join(os.path.expanduser("~"), "data/Aquarius/output")
ro = RamsesOutput(datadir, ioutput)


# Map operator : minimum temperature along line-of-sight
class MyTempOperator(Operator):
    def __init__(self):
        def invT_func(dset):
            P = dset["P"]
            rho = dset["rho"]
            r = rho / P
            # print r[(rho<=0.0)+(P<=0.0)]
            # r[(rho<=0.0)*(P<=0.0)] = 0.0
            return r

        d = {"invTemp": invT_func}
        Operator.__init__(self, d, is_max_alos=True)

    def operation(self, int_dict):
        map = list(int_dict.values())[0]
        mask = (map == 0.0)
        mask2 = map != 0.0
        map[mask2] = 1.0 / map[mask2]
        map[mask] = 0.0
        return map


scal_op = MyTempOperator()

# Map region
center = [0.567111, 0.586555, 0.559156]
axes = {"los": "z"}

# Map processing
rt = raytracing.RayTracer(ro, ["rho", "P"])
for axname, axis in axes.items():
    cam = Camera(center=center, line_of_sight_axis=axis, up_vector="y", region_size=[3.0E-3, 3.0E-3],
                 size_unit=ro.info["unit_length"], distance=1.5E-3, far_cut_depth=1.5E-3, map_max_size=512)
    map = rt.process(scal_op, cam)
    factor = ro.info["unit_temperature"].express(C.K)
    scale = ro.info["unit_length"].express(C.Mpc)

    # Save map into HDF5 file
    mapname = "gas_rt_Tmin_%s_%5.5i" % (axname, ioutput)
    h5fname = ImgPlot.save_map_HDF5(map, cam, map_name=mapname)

    # Plot map into Matplotlib figure/PIL Image
    fig = ImgPlot.save_HDF5_to_plot(h5fname, map_unit=("K", factor), axis_unit=("Mpc", scale), cmap="hot", fraction=0.0)
    # pil_img = save_HDF5_to_img(h5fname, cmap="hot")

# Save into PNG image file
# save_HDF5_to_plot(h5fname, map_unit=("K",factor), axis_unit=("Mpc", scale), img_path="./", cmap="hot")
# save_HDF5_to_img(h5fname, img_path="./", cmap="hot")

# pylab.show()
