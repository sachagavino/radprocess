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

from pymses.utils import constants as C
from nose import tools
import os


class TestConstantsUnit:
    _h5f = None
    _h5f_v1 = None

    @classmethod
    def setup_class(cls):
        import tables as T
        import h5py
        PWD = os.path.dirname(__file__)
        old_dir = os.path.join(PWD, "pymses4_old")
        units_v1_dir = os.path.join(PWD, "units_v1")
        cls._h5f = T.open_file(os.path.join(old_dir, "units.h5"), 'r')
        cls._h5f_v1 = h5py.File(os.path.join(units_v1_dir, "units.h5"), 'r')

    def test_old_surface_density_unit(self):
        print("Test old format surface density unit")
        surf_density_group = self._h5f.get_node("/unit_surf_density")
        sigma_unit = C.Unit._load_legacy_HDF5_unit(surf_density_group)
        tools.assert_equal(sigma_unit, 10.0 * C.Msun / C.pc ** 2)

    def test_old_velocity_unit(self):
        print("Test old format velocity unit")
        unit_velocity_group = self._h5f.get_node("/unit_velocity")
        vel_unit = C.Unit._load_legacy_HDF5_unit(unit_velocity_group)
        tools.assert_equal(vel_unit, 65.0 * C.km / C.s)

    def test_old_magnetic_field_unit(self):
        print("Test old format magnetic field unit")
        unit_magn_group = self._h5f.get_node("/unit_magn")
        B_unit = C.Unit._load_legacy_HDF5_unit(unit_magn_group)
        tools.assert_equal(B_unit, 4.2 * C.mGauss)

    def test_old_temperature_unit(self):
        print("Test old format temperature unit")
        unit_temp_group = self._h5f.get_node("/unit_temp")
        temp_unit = C.Unit._load_legacy_HDF5_unit(unit_temp_group)
        tools.assert_equal(temp_unit, 273.15 * C.K)

    def test_luminosity_unit_v1(self):
        print("Test v1 luminosity unit")
        lum_unit = C.Unit.from_HDF5(self._h5f_v1['/catalog/unit_luminosity'])
        tools.assert_equal(lum_unit, 3.5 * C.Lsun)
        tools.assert_equal(lum_unit.name, "3.5 Lsun")
        tools.assert_equal(lum_unit.description, "Reference luminosity : 3.5 solar luminosity")

    def test_force_unit_v1(self):
        print("Test v1 force unit")
        force_unit = C.Unit.from_HDF5(self._h5f_v1['/catalog/unit_force'])
        tools.assert_equal(force_unit, 53.1 * C.N)
        tools.assert_equal(force_unit.name, "53.1 N")
        tools.assert_equal(force_unit.description, "Reference force unit : 53.1 Newton")

    def test_energy_unit_v1(self):
        print("Test v1 energy unit")
        engy_unit = C.Unit.from_HDF5(self._h5f_v1['/catalog/unit_energy'])
        tools.assert_equal(engy_unit, 5.0 * C.erg)
        tools.assert_equal(engy_unit.name, "5 erg")
        tools.assert_equal(engy_unit.description, "Reference energy value : 5 erg")

    def test_angle_unit_v1(self):
        print("Test v1 angle unit")
        angle_unit = C.Unit.from_HDF5(self._h5f_v1['/catalog/unit_angle'])
        tools.assert_equal(angle_unit, 2 * C.hourangle)
        tools.assert_equal(angle_unit.name, "2 hourangle")
        tools.assert_equal(angle_unit.description, "Reference angle : pi / 6")

    def test_micron_unit_v1(self):
        print("Test v1 micron unit")
        um_unit = C.Unit.from_HDF5(self._h5f_v1['/catalog/micron_unit'])
        tools.assert_equal(um_unit, C.um)
        tools.assert_equal(um_unit.name, "um")
        tools.assert_equal(um_unit.description, "Micron")
        tools.assert_equal(um_unit.latex, "$\mu m$")

    @classmethod
    def teardown_class(cls):
        if cls._h5f is not None:
            cls._h5f.close()
        if cls._h5f_v1 is not None:
            cls._h5f_v1.close()
