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
:mod:`pymses.utils.constants.units` --- physical unit module
============================================================

"""
from numpy import isscalar, array, abs, isclose, log10, argmin
from .physical_type import get_phystype
from pymses.utils import HDF5Serializable


class UnitError(Exception):
    pass


_unit_registry = {}


class Unit(HDF5Serializable):
    r"""
    Dimensional physical unit class

    Parameters
    ----------
    name : ``string``
        Unit name
    base_unit : ``Unit`` instance
        Composite unit from which this instance should be initialised
    coeff  : ``float``
        dimensionless value of the unit instance.
    dims : 8-``tuple`` of ``int``
        dimension of the unit object expressed in the international unit system
        (``m``, ``kg``, ``s``, ``K``, ``rad``, ``A``, ``mol``, ``cd``)
    descr : ``string`` or `None`
        Unit description
    latex: ``string`` or `None`
        Unit displayed name (latex format)

    Examples
    --------

        >>> cs_m_s = Unit(name="m/s", coeff=340.0, dims=(1, 0, -1, 0, 0, 0, 0, 0), descr="velocity SI unit")
        >>> print "sound speed = %g m/h"%(cs_m_s.express(km/hour))
        sound speed = 1224 km/h
        >>>
        >>> dens = Unit(name="%s.%s$^{-3}$" % (Msun.latex, kpc.latex), base_unit=C.Msun/C.kpc**3,
        ... descr="Solar mass per cubic kiloparsec")
        >>> print dens
        (6.76957356533e-29 m^-3.kg)
    """
    _unit_system = ["m", "kg", "s", "K", "rad", "A", "mol", "cd"]

    def __new__(cls, name="", base_unit=None, coeff=1.0, dims=None, descr=None, latex=None):
        if name != "" and name in _unit_registry:
            return _unit_registry[name]

        return super(Unit, cls).__new__(cls)

    def __init__(self, name="", base_unit=None, coeff=1.0, dims=None, descr=None, latex=None):
        # Already initialised Unit instance (in the registry)
        if name != "" and name in _unit_registry:
            return

        # Unit name + laTex displayed name
        self._name = name
        self._latex = latex

        if base_unit is not None:
            if not isinstance(base_unit, Unit):
                raise AttributeError("'base_unit' parameter must be a a 'Unit' instance.")
            self._dimensions = array(base_unit._dimensions, 'i')
            self._coeff = base_unit._coeff
        elif dims is None:
            raise AttributeError("'dims' attribute is mandatory to define a 'Unit' object.")
        else:
            if len(dims) != len(Unit._unit_system):
                raise AttributeError("Wrong dimensions : must be a tuple of length " + str(len(Unit._unit_system)))
            self._dimensions = array(dims, 'i')
            self._coeff = coeff

        if descr is None:
            self._descr = ""
        else:
            self._descr = descr

    @classmethod
    def from_name(cls, unit_name):
        if unit_name not in _unit_registry:
            raise AttributeError("Unknown unit name '%s'" % unit_name)

        return _unit_registry[unit_name]

    def __eq__(self, other):
        """

        Parameters
        ----------
        other: :class:`~pymses.utils.constants.unit.Unit`
            other Unit instance to compare to

        Returns
        -------
        e: ``bool``
            True if coeff and dimenions are identical, otherwise False.
        """
        if not isclose(self._coeff, other._coeff, rtol=1.0e-6):
            return False
        if ((self._dimensions - other.dimensions) != 0).any():
            return False
        return True

    @property
    def name(self):
        """
        Unit name
        """
        return self._name

    @property
    def latex(self):
        """
        Unit name (LaTex format)
        """
        if self._latex is not None:
            return self._latex
        return self._name

    @property
    def dimensions(self):
        """
        Unit dimension array
        """
        return self._dimensions

    @property
    def description(self):
        """
        Unit description
        """
        return self._descr

    @property
    def coeff(self):
        """
        Constant value of this unit
        :return:
        """
        return self._coeff

    def is_base_unit(self):
        if self._coeff != 1.0:
            return False

        msk = (self._dimensions != 0)
        n = msk.sum()
        if n == 0:  # Dimensionless unit
            return True
        elif n == 1:
            if self._dimensions[msk][0] == 1:
                return True

        return False

    def info(self):
        """
        Print information about this Unit object.

        If any, print the name and description of this unit, then print the value of this unit and the list
        of equivalent unit contained in the built-in unit registry associated with their conversion factor.
        """
        s = ""
        if self._name != "":
            title = "Unit : %s" % self._name
            l = len(title)
            s += title + "\n" + (l * "-") + "\n"

        # Description
        if len(self._descr) > 0:
            s += self._descr + "\n\n"

        s += "Value\n" \
             "-----\n" + str(self._coeff) + " " + self._decompose_base_units() + "\n\n"

        s += "Equivalent units\n" \
             "----------------\n" + ""

        for u in self.equivalent_unit_list():
            s += " * %-12s : %14g %s\n" % (u.name, u.express(self), self._name)

        print(s)

    def equivalent_unit_list(self):
        """
        Get the equivalent unit list (with same physical type)
        :return: equivalent Unit object list
        """
        l = []
        t = self.physical_type
        for u in _unit_registry.values():
            if self == u:
                continue
            if u.physical_type == t:
                l.append(u)

        return l

    def appropriate_unit(self, nearest_log10=1.0):
        """
        Try to find the better suited unit (among available equivalent units to represent this unit
        """
        eq_units = self.equivalent_unit_list()
        if len(eq_units) == 0:
            return self

        tlist = [(eq_units, self.express(eq_unit)) for eq_unit in eq_units]
        imin = argmin(array([abs(log10(t[1]) - nearest_log10) for t in tlist]))
        best_unit, ph_val = tlist[imin]
        return ph_val, best_unit

    @property
    def physical_type(self):
        """
        Get the unit physical type (dimensioned physical quantity).
        :return:
        """
        return get_phystype(self)

    def __div__(self, other):
        r"""
        Returns the result of the division of `self` by `other`

        """
        if isscalar(other):
            newdims = self._dimensions
            newval = self._coeff / other
            return Unit(coeff=newval, dims=newdims)
        elif isinstance(other, Unit):
            newdims = self._dimensions - other._dimensions
            newval = self._coeff / other._coeff
            return Unit(coeff=newval, dims=newdims)
        else:
            raise UnitError("Unable to divide a Unit instance by something which is neither a Unit object nor a "
                            "scalar.")

    def __truediv__(self, other):
        r"""
        Returns the result of the future division of `self` by `other`

        """
        return self.__div__(other)

    def __rdiv__(self, other):
        r"""
        Returns the result of the division of `other` by `self`

        """
        if isscalar(other):
            newdims = -self._dimensions
            newval = other / self._coeff
            return Unit(coeff=newval, dims=newdims)
        elif isinstance(other, Unit):
            newdims = other._dimensions - self._dimensions
            newval = other._coeff / self._coeff
            return Unit(coeff=newval, dims=newdims)
        else:
            raise UnitError("Unable to divide something which is neither a Unit object nor a scalar by a Unit "
                            "instance.")

    def __rtruediv__(self, other):
        r"""
        Returns the result of the future division of `other` by `self`

        """
        return self.__rdiv__(other)

    def __mul__(self, other):
        r"""
        Returns the result of the multiplication of `self` with `other`

        """
        if isscalar(other):
            newdims = self._dimensions
            newval = self._coeff * other
            return Unit(coeff=newval, dims=newdims)
        elif isinstance(other, Unit):
            newdims = self._dimensions + other._dimensions
            newval = self._coeff * other._coeff
            return Unit(coeff=newval, dims=newdims)
        else:
            raise UnitError("Unable to multiply a Unit instance by something which is neither a Unit object nor a "
                            "scalar.")

    def __rmul__(self, other):
        r"""
        Returns the result of the multiplication of `other` with `self`

        """
        return self * other

    def __pow__(self, n):
        r"""
        Returns the result of the exponentiation of `self` by `n`

        """
        return Unit(coeff=self._coeff ** n, dims=n * self._dimensions)

    def _decompose_base_units(self):
        s = ""

        d = array(self._dimensions)
        for i in range(len(Unit._unit_system)):
            if d[i] != 0:
                if (d[0:i] != 0).any():
                    s += "."
                s += Unit._unit_system[i]
                if d[i] != 1:
                    s += "^" + str(d[i])
        return s

    def __repr__(self):
        r"""
        ``string`` representation of `self`

        """
        strrep = ""

        # Name & value
        if len(self._name) > 0:
            strrep += self._name + " : "

        strrep += "(" + str(self._coeff) + " " + self._decompose_base_units() + ")"
        return strrep

    def express(self, unit):
        r"""
        Unit conversion method. Gives the conversion factor of this ``Unit`` object
        expressed into another (dimension-compatible) given `unit`.

        Checks that :

        * the `unit` param. is also a ``Unit`` instance
        * the `unit` param. is dimension-compatible

        Parameters
        ----------
        unit : ``Unit`` object
            unit in which the conversion is made

        Returns
        -------
        fact : ``float``
            conversion factor of `self` expressed in `unit`

        Examples
        --------
        * Conversion of a kpc expressed in light-years :

            >>> factor = kpc.express(ly)
            >>> print "1 kpc = %f ly"%factor
            1 kpc = 3261.563163 ly

        * Conversion of :math:`1 M_{\odot}` into kpc/Myr :

            >>> print Msun.express(kpc/Myr)
            Incompatible dimensions between :
            - Msun : (1.9889e+30 kg) (type: mass) and
            - (977792.221673 m.s^-1) (type: velocity)


        """
        if not isinstance(unit, Unit):
            raise AttributeError("'unit' attribute must be a 'Unit' instance.")

        dmax = len(Unit._unit_system) - 1
        if tuple(self._dimensions[0:dmax]) != tuple(unit._dimensions[0:dmax]):
            raise UnitError("Incompatible dimensions between :\n"
                            "- " + str(self) + " (type: " + self.physical_type + ") and\n" +
                            "- " + str(unit) + " (type: " + unit.physical_type + ")")
        return (self / unit)._coeff

    def _h5_serialize(self, h5group, **kwargs):
        """
        Serialize the Unit instance into a HDF5 file.

        Parameters
        ----------
        h5group: ``h5py.highlevel.Group``
            Main group to serialize the Unit into.
        """
        h5group.attrs['name'] = str(self._name)
        h5group.create_dataset("dimensions", data=self._dimensions)
        h5group.attrs['coefficient'] = self._coeff
        h5group.attrs['description'] = str(self._descr)
        if self._latex is not None:
            h5group.attrs['latex'] = str(self._latex)

    @classmethod
    def _h5_deserialize(cls, h5group, version):
        """
        Read  a unit from a HDF5 file.

        Parameters
        ----------
        h5group: ``h5py.highlevel.Group``
            Main group to deserialize the Unit from.
        version: ``int``
            Version of the Unit class to deserialize.

        Returns
        -------
        cam: :class:`~pymses.utils.constants.unit.Unit`
            new Unit instance
        """
        if version == Unit._pymses_h5_version:  # Current Unit version
            name = h5group.attrs['name']
            udims = h5group['dimensions'][:]
            ucoeff = h5group.attrs['coefficient']
            descr = h5group.attrs['description']
            lat = None
            if "latex" in h5group.attrs:
                lat = h5group.attrs["latex"]
            return Unit(name, coeff=ucoeff, dims=udims, descr=descr, latex=lat)
        else:
            raise AttributeError("Unknown Unit version (%d)" % int(version))

    @classmethod
    def _load_legacy_HDF5_unit(cls, h5g):
        """
        Load a legacy (Pymses v<=4.0.0) unit object from HDF5 group

        Parameters
        ----------
        h5g: :class:`tables.H5Group`

        Returns
        -------
        u: :class:`~pymses.utils.constants.unit.Unit`
            Unit object
        """
        # Old dimensions are in (m, kg, s, K, h, T) while current version Unit dimensions are in
        # (m, kg, s, K, rad, A, mol, cd) => the Tesla unit must be converted and the h unit need to be abandoned
        olddims = h5g.dimensions.read()
        Tpow = olddims[5]
        # 1 Tesla = 1 kg/(s**2 * A)
        udims = array([olddims[0], olddims[1] + Tpow, olddims[2] - 2 * Tpow, olddims[3], 0, -Tpow, 0, 0], dtype='i')
        uval = h5g.val.read()
        return Unit(coeff=uval, dims=udims)


def create_unit(name="", base_unit=None, coeff=1.0, dims=None, descr=None, latex=None):
    """
    Add a new Unit instance to the registry
    """
    if name in _unit_registry:
        raise ValueError("%s Unit is already defined as %s" % (name, str(_unit_registry[name])))

    u = Unit(name, base_unit, coeff, dims, descr, latex)
    _unit_registry[name] = u
    return u


__all__ = ["Unit"]
