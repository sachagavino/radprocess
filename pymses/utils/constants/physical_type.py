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
:mod:`pymses.utils.constants.physical_type` --- Physical quantity type module.
==============================================================================

"""

_phystype_unit_mapping = {}


def create_phystype(unit, name):
    """
    Adds a new physical unit mapping.

    Parameters
    ----------
    unit : `~astropy.units.UnitBase` instance
        The unit to map from.

    name : str
        The physical name of the unit.
    """
    d = tuple(unit.dimensions)
    if d in _phystype_unit_mapping:
        raise ValueError("%s (%s) already defined as %s" % (str(d), name, _phystype_unit_mapping[d]))
    _phystype_unit_mapping[d] = name


def get_phystype(unit):
    """
    Given a unit, returns the name of the physical quantity it represents.  If it represents an unknown physical
    quantity, 'unknown' is returned.

    Parameters
    ----------
    unit : ``Unit`` instance
        The requiered unit object

    Returns
    -------
    t : ``string`
        The name of the physical quantity, or 'unknown' if not known.
    """
    d = tuple(unit.dimensions)
    return _phystype_unit_mapping.get(d, 'unknown')

__all__ = ['get_phystype']
