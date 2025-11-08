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


import numpy as N
from PIL import Image as I
import os
from pymses.analysis import DataMap
from pymses.analysis.plot import PlainPNGImage
from pymses.utils import constants as C
from nose import tools
import tempfile


class TestDatamapHDF5:
    _h5f = None

    @classmethod

    def setup_class(cls):
        import h5py
        cls._pwd = os.path.dirname(__file__)
        cls._h5_dir = os.path.join(cls._pwd, "pymses4_old")
        cls._h5_dir_v1 = os.path.join(cls._pwd, "datamap_v1")
        cls._h5_dir_v2 = os.path.join(cls._pwd, "datamap_v2")
        cls._h5fname_legacy = os.path.join(cls._h5_dir, "Aquarius_density_RT_3disks.h5")
        cls._h5f_v1 = h5py.File(os.path.join(cls._h5_dir_v1, "Aquarius_density_RT_3disks_v1.h5"), 'r')
        cls._h5f_v2 = h5py.File(os.path.join(cls._h5_dir_v2, "Aquarius_density_RT_3disks_v2.h5"), 'r')
        cls._datamap_png_image_v2 = os.path.join(cls._h5_dir_v2, "Aquarius_density_RT_3disks.png")
        fd, cls._h5f_v2_temp = tempfile.mkstemp(prefix="pymses_datamap", suffix=".h5")
        os.close(fd)
        fd, cls._temp_img = tempfile.mkstemp(prefix="pymses_img", suffix=".png")
        os.close(fd)

    def test_datamap_legacy_to_v1_hdf5_io(self):
        print("Test Datamap object HDF5 I/O backward compatibility (legacy => v1)")
        d_legacy = DataMap.load_legacy_HDF5_datamap(self._h5fname_legacy)
        d1 = DataMap.from_HDF5(self._h5f_v1)
        tools.assert_equal(d_legacy, d1)

    def test_datamap_v1_to_2_hdf5_io(self):
        print("Test Datamap object HDF5 I/O backward compatibility (v1 => v2)")
        d1 = DataMap.from_HDF5(self._h5f_v1)
        d1.save_HDF5(self._h5f_v2_temp, float32=True)

        # Density custom unit
        # mu = C.Unit(name="%s/%s**3" % (C.Msun.name, C.kpc.name), base_unit=C.Msun/C.kpc**3,
        #             descr="Solar mass per cubic kiloparsec", latex="%s.%s$^{-3}$" % (C.Msun.latex, C.kpc.latex))

        d2_temp = DataMap.from_HDF5(self._h5f_v2_temp)
        tools.assert_equal(d1, d2_temp)

        # img_fname = os.path.join(self._h5_dir_v1, "Aquarius_density_RT_3disks.png")
        # plot_fname = os.path.join(self._h5_dir_v1, "Aquarius_density_RT_3disks_plot.png")
        # d2.save_PNG(img_fname, fraction=10.0, cmap="BlackPurpRedBlGreen", verbose=True)
        # d2.save_plot(plot_fname, fraction=10.0, cmap="BlackPurpRedBlGreen", axis_unit=C.Mpc, map_unit=mu, verbose=True)

        d2 = DataMap.from_HDF5(self._h5f_v2)
        tools.assert_equal(d2, d2_temp)

    def test_datamap_legacy_to_png_image(self):
        print("Test Datamap object HDF5 I/O backward compatibility (legacy => PNG image)")
        d_legacy = DataMap.load_legacy_HDF5_datamap(self._h5fname_legacy)

        img_generator = PlainPNGImage(cmap="BlackPurpWhiBlGreen")
        img = d_legacy.save_PNG(img_generator, img_fname=self._temp_img, fraction=10.0, verbose=True)

        img2 = I.open(self._datamap_png_image_v2, 'r')
        aimg = N.array(img.getdata())
        aimg2 = N.array(img2.getdata())
        img2.close()
        tools.assert_true(N.allclose(aimg, aimg2, rtol=1.0e-6))

    def test_datamap_v1_to_png_image(self):
        print("Test Datamap object HDF5 I/O backward compatibility (v1 => PNG image)")
        d1 = DataMap.from_HDF5(self._h5f_v1)

        img_generator = PlainPNGImage(cmap="BlackPurpWhiBlGreen")
        img = d1.save_PNG(img_generator, img_fname=self._temp_img, fraction=10.0, verbose=True)

        img2 = I.open(self._datamap_png_image_v2, 'r')
        aimg = N.array(img.getdata())
        aimg2 = N.array(img2.getdata())
        img2.close()
        tools.assert_true(N.allclose(aimg, aimg2, rtol=1.0e-6))

    def test_datamap_v2_to_png_image(self):
        print("Test Datamap object HDF5 I/O backward compatibility (v2 => PNG image)")
        d2 = DataMap.from_HDF5(self._h5f_v2)

        img_generator = PlainPNGImage(cmap="BlackPurpWhiBlGreen")
        img = d2.save_PNG(img_generator, img_fname=self._temp_img, fraction=10.0, verbose=True)

        img2 = I.open(self._datamap_png_image_v2, 'r')
        aimg = N.array(img.getdata())
        aimg2 = N.array(img2.getdata())
        img2.close()
        tools.assert_true(N.allclose(aimg, aimg2, rtol=1.0e-6))

    @classmethod
    def teardown_class(cls):
        if cls._h5f_v1 is not None:
            cls._h5f_v1.close()
        if cls._h5f_v2 is not None:
            cls._h5f_v2.close()
        if os.path.isfile(cls._temp_img):
            os.remove(cls._temp_img)
        if os.path.isfile(cls._h5f_v2_temp):
            os.remove(cls._h5f_v2_temp)
