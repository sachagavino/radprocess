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
:mod:`pymses.analysis.quantities.generic` --- Generic quantity module
---------------------------------------------------------------------

"""
from .base import ScalarQuantity, VectorQuantity, MultivaluedQuantity
from pymses.utils import constants as C


class DensityQty(ScalarQuantity):
    """
    Density scalar quantity

    Parameters
    ----------
    name: ``string``
        Scalar quantity name. Default "rho".
    amr_field_list: ``list``or None.
        AMR field name list. Default None : use ['rho'].
    tex_name: ``string``
        LaTex display quantity name. Default r"\rho".
    """
    def __init__(self, name="rho", amr_field_list=None, unit=None, tex_name=r"\rho"):
        if amr_field_list is None:
            amr_fields = ["rho"]
        else:
            amr_fields = amr_field_list
        if unit is None:
            u = C.g_cc
        else:
            u = unit
        super(DensityQty, self).__init__(name, amr_fields, u, tex_name)

    def func_dset(self, amr_dataset):
        return amr_dataset["rho"]

    def data_unit(self, info):
        return info["unit_density"]


class TemperatureQty(ScalarQuantity):
    """
    Temperature scalar quantity

    Parameters
    ----------
    name: ``string``
        Scalar quantity name. Default "T".
    amr_field_list: ``list``or None.
        AMR field name list. Default None : use ['P', 'rho'].
    tex_name: ``string``
        LaTex display quantity name. Default r"\textrm{T}".
    """
    def __init__(self, name="T", amr_field_list=None, unit=None, tex_name=r"\textrm{T}"):
        if amr_field_list is None:
            amr_fields = ["rho", "P"]
        else:
            amr_fields = amr_field_list
        if unit is None:
            u = C.K
        else:
            u = unit
        super(TemperatureQty, self).__init__(name, amr_fields, u, tex_name)

    def func_dset(self, amr_dataset):
        return amr_dataset["P"]/amr_dataset["rho"]

    def data_unit(self, info):
            return info["unit_temperature"] * info["mu_gas"]


class MassQty(ScalarQuantity):
    """
    Mass scalar quantity

    Parameters
    ----------
    name: ``string``
        Scalar quantity name. Default "mass".
    amr_field_list: ``list``or None.
        AMR field name list. Default None : use ['rho'].
    tex_name: ``string``
        LaTex display quantity name. Default r"\textrm{mass}".
    """
    def __init__(self, name="mass", amr_field_list=None, unit=None, tex_name=r"\textrm{mass}"):
        if amr_field_list is None:
            amr_fields = ["rho"]
        else:
            amr_fields = amr_field_list
        if unit is None:
            u = C.Msun
        else:
            u = unit
        super(MassQty, self).__init__(name, amr_fields, u, tex_name)

    def func_dset(self, amr_dataset):
        if 'size' in amr_dataset.fields:
            mass = amr_dataset["rho"] * amr_dataset["size"]**3
        elif 'level' in amr_dataset.fields:
            mass = amr_dataset["rho"] * 1./2**(3*amr_dataset["level"])
        else:
            raise ValueError("Cannot compute AMR cell size from this dataset: 'size' or 'level' field is missing.")
        return mass

    def data_unit(self, info):
            return info["unit_mass"]


class PhiQty(ScalarQuantity):
    """
    Gravitational potential scalar quantity

    Parameters
    ----------
    name: ``string``
        Scalar quantity name. Default "phi".
    amr_field_list: ``list``or None.
        AMR field name list. Default None : use ['phi'].
    tex_name: ``string``
        LaTex display quantity name. Default r"\phi".
    """
    def __init__(self, name="phi", amr_field_list=None, unit=None, tex_name=r"\phi"):
        if amr_field_list is None:
            amr_fields = ["phi"]
        else:
            amr_fields = amr_field_list
        if unit is None:
            u = C.Unit(base_unit=(C.m/C.s)**2, descr="Squared meters per square second", latex=r"m$^{2}$.s$^{-2}$")
        else:
            u = unit
        super(PhiQty, self).__init__(name, amr_fields, u, tex_name)

    def func_dset(self, amr_dataset):
        return amr_dataset["phi"]

    def data_unit(self, info):
            return info["unit_gravpot"]


class VelocityQty(VectorQuantity):
    """
    Velocity vector quantity

    Parameters
    ----------
    name: ``string``
        Vector quantity name. Default "V".
    amr_field_list: ``list``or None.
        AMR field name list. Default None : use ['vel'].
    tex_name: ``string``
        LaTex display quantity name. Default r"\textrm{v}".
    """
    def __init__(self, name="vel", amr_field_list=None, unit=None, tex_name=r"\textrm{v}"):
        if amr_field_list is None:
            amr_fields = ["vel"]
        else:
            amr_fields = amr_field_list
        if unit is None:
            u = C.km_s
        else:
            u = unit
        super(VelocityQty, self).__init__(name, amr_fields, u, tex_name)

    def func_dset(self, amr_dataset):
        return amr_dataset["vel"]

    def data_unit(self, info):
        return info["unit_velocity"]


class BQty(VectorQuantity):
    """
    Magnetic field vector quantity

    Parameters
    ----------
    name: ``string``
        Vector quantity name. Default "B".
    amr_field_list: ``list``or None.
        AMR field name list. Default None : use ['Bl', 'Br'].
    tex_name: ``string``
        LaTex display quantity name. Default r"\textrm{B}".
    """
    def __init__(self, name="B", amr_field_list=None, unit=None, tex_name=r"\textrm{B}"):
        if amr_field_list is None:
            amr_fields = ["vel"]
        else:
            amr_fields = amr_field_list
        if unit is None:
            u = C.Gauss
        else:
            u = unit
        super(BQty, self).__init__(name, amr_fields, u, tex_name)

    def func_dset(self, amr_dataset):
        return (amr_dataset['Bl'] + amr_dataset['Br']) / 2.

    def data_unit(self, info):
        return info['unit_mag']


class GravAccQty(VectorQuantity):
    """
    Gravitational acceleration field vector quantity

    Parameters
    ----------
    name: ``string``
        Vector quantity name. Default "g".
    amr_field_list: ``list``or None.
        AMR field name list. Default None : use ['g'].
    tex_name: ``string``
        LaTex display quantity name. Default r"\textrm{g}".
    """
    def __init__(self, name="g", amr_field_list=None, unit=None, tex_name=r"\textrm{g}"):
        if amr_field_list is None:
            amr_fields = ["g"]
        else:
            amr_fields = amr_field_list
        if unit is None:
            u = C.Unit(base_unit=C.m/C.s**2, descr="Meters per square second", latex=r"m.s$^{-2}$")
        else:
            u = unit
        super(GravAccQty, self).__init__(name, amr_fields, u, tex_name)

    def func_dset(self, amr_dataset):
        return amr_dataset['g']

    def data_unit(self, info):
        return info["unit_velocity"] / info["unit_time"]  # * info_dict["boxlen"] ?


class ErQty(MultivaluedQuantity):
    """
    Multigroup energy field multivalued quantity

    Parameters
    ----------
    name: ``string``
        Multivalued quantity name. Default "Er".
    amr_field_list: ``list``or None.
        AMR field name list. Default None : use ['Er'].
    tex_name: ``string``
        LaTex display quantity name. Default r"\textrm{E_{r}}".
    """
    def __init__(self, name="Er", amr_field_list=None, unit=None, tex_name=r"\textrm{E_{r}}"):
        if amr_field_list is None:
            amr_fields = ["Er"]
        else:
            amr_fields = amr_field_list
        if unit is None:
            u = C.Unit(base_unit=C.erg/C.cm**3, descr="Ergs per cubic centimeters",
                       latex=r"$\textrm{erg}.\textrm{cm}^{-3}$")
        else:
            u = unit
        super(ErQty, self).__init__(name, amr_fields, u, tex_name)

    def func_dset(self, amr_dataset):
        return amr_dataset['Er']

    def data_unit(self, info):
        return info["unit_density"] * info["unit_velocity"]**2


__all__ = ['DensityQty', "TemperatureQty", 'MassQty', 'VelocityQty', 'GravAccQty', 'PhiQty', 'ErQty', 'BQty']
