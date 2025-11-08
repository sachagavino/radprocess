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
:mod:`pymses.analysis.operator` --- post-processing operator package
--------------------------------------------------------------------

"""
from numpy import zeros_like, arctan2

from pymses.core import IsotropicExtPointDataset
from pymses.utils import constants as C


class AbstractOperator(object):
    r"""
    Abstract operator generic class

    Parameters
    ----------
    scalar_field_dict: ``string``-indexed ``dict`` of ``callable``
        scalar field function dictionary
    output_unit: :class:`~pymses.utils.constants.unit.Unit`
        output data unit
    is_max_alos: ``boolean``
        if True, take the maximum along the line-of-sight of each scalar field. Otherwise, take the sum of the scalar
        fields integrated along the line-of-sight. Default False.
    use_cell_dx: ``boolean``
        Need to use the cell size ? Default False.
    """

    def __init__(self, scalar_field_dict, output_unit, is_max_alos=False, use_cell_dx=False):
        if not isinstance(scalar_field_dict, dict):
            raise AttributeError("'scalar_func_dict' must be a dictionary.")

        if not isinstance(output_unit, C.unit.Unit):
            raise AttributeError(output_unit + " is not a valid Unit instance")

        self._scalar_func_dict = {}
        self._output_unit = output_unit

        for key, t in scalar_field_dict.items():
            if not isinstance(key, str) or not hasattr(t, '__call__'):
                raise AttributeError("'scalar_field_dict' must contain string-indexed callables")

        self._scalar_func_dict = scalar_field_dict
        self._max_alos = is_max_alos
        self._use_cell_size = use_cell_dx

    def __iter__(self):
        return self.iter_scalar_func()

    def iter_scalar_func(self):
        """
        TODO

        Returns
        -------
        """
        for key, func in self._scalar_func_dict.items():
            yield key, func

    def varnames(self):
        """
        TODO

        Returns
        -------
        """
        return list(self._scalar_func_dict.keys())

    @property
    def output_unit(self):
        return self._output_unit

    def nscal_func(self):
        """
        Returns the total number of scalar field

        Returns
        -------
        nscal_func: ``int``
            total number of scalar fields.
        """
        return len(self._scalar_func_dict)

    def is_max_alos(self):
        """
        TODO

        Returns
        -------
        """
        return self._max_alos

    def use_cell_dx(self):
        """
        TODO

        Returns
        -------
        """
        return self._use_cell_size

    def operation(self, maps):
        """
        TODO

        Parameters
        ----------
        maps:

        Returns
        -------

        """
        raise NotImplementedError()


class MultiFieldOperator(AbstractOperator):
    """
    TODO
    """

    def __init__(self, scalar_field_dict, output_unit_dict, **kwargs):
        first_unit = None
        if not isinstance(output_unit_dict, dict) or len(output_unit_dict) == 0:
            raise AttributeError(output_unit_dict + " is not a valid Unit dictionary")
        else:
            for unit in output_unit_dict.values():
                if not isinstance(unit, C.unit.Unit):
                    raise AttributeError(unit + " is not a valid unit")
                if first_unit is None:
                    first_unit = unit

        super(MultiFieldOperator, self).__init__(scalar_field_dict, first_unit, **kwargs)
        # Override private _output_unit attribute by the unit dictionary
        self._output_unit = output_unit_dict

    def operation(self, maps):
        raise NotImplementedError()


class ScalarOperator(AbstractOperator):
    r"""
    ScalarOperator class

    Parameters
    ----------
    scalar_func : ``function``
        single `dset` argument function returning the scalar data ``array`` from this `dset` AbstractDataset.
    scalar_unit : :class:`~pymses.utils.constants.unit.Unit`
        scalar field unit.

    Examples
    --------
    >>> # Density field scalar operator
    >>> op = ScalarOperator(lambda dset: dset["rho"], ro.info["unit_density"])

    """

    def __init__(self, scalar_func, scalar_unit, **kwargs):
        super(ScalarOperator, self).__init__({"scalar": scalar_func}, scalar_unit, **kwargs)

    def operation(self, maps):
        return maps["scalar"]


class FractionOperator(AbstractOperator):
    r"""
    FractionOperator class

    Parameters
    ----------
    up_func   : ``function``
        numerator function like `scalar_func` (see :class:`~pymses.analysis.ScalarOperator`)
    down_func : ``function``
        denominator function like `scalar_func` (see :class:`~pymses.analysis.ScalarOperator`)

    Examples
    --------
    >>> # Mass-weighted pressure scalar operator
    >>> num = lambda dset: dset["rho"]*dset["P"]   * dset.get_sizes()**3
    >>> den = lambda dset: dset["rho"] * dset.get_sizes()**3
    >>> op = FractionOperator(num, den)

    .. math::
        I = \frac{\int\limits_{V} \rho \times P \mathrm{d}V}{\int\limits_{V} \rho \mathrm{d}V}

    """

    def __init__(self, num_func, denom_func, frac_unit, **kwargs):
        super(FractionOperator, self).__init__({"numerator": num_func, "denominator": denom_func}, frac_unit, **kwargs)

    def operation(self, maps):
        md = maps["denominator"]
        mu = maps["numerator"]
        map = zeros_like(md)
        mask = (md != 0.0)
        map[mask] = mu[mask] / md[mask]
        return map


class PerpendicularVectorFieldDirectionOperator(AbstractOperator):
    r"""
    PerpendicularVectorFieldDirectionOperator class

    Parameters
    ----------
    vec_x_func   : ``function``
        vector component x-axis component function like `scalar_func` (see :class:`~pymses.analysis.ScalarOperator`)
    vec_y_func : ``function``
        vector component y-axis component function like `scalar_func` (see :class:`~pymses.analysis.ScalarOperator`)

    Examples
    --------
    >>> # Mass-weighted velocity direction
    >>> camaxes = camera.get_camera_axis()
    >>> vel_x = lambda dset: dset["rho"] * numpy.dot(dset["vel"], camaxes[0, :])
    >>> vel_y = lambda dset: dset["rho"] * numpy.dot(dset["vel"], camaxes[1, :])
    >>> velocity_dir_xy_plane_op = PerpendicularVectorFieldDirectionOperator(vel_x, vel_y)

    .. math::
        Ix = \int\limits_{V} \rho \times Vx \mathrm{d}V}
        Iy = \int\limits_{V} \rho \times Vy \mathrm{d}V}
        alpha = arctan2(Iy, Ix)

    """

    def __init__(self, vec_x_func, vec_y_func, **kwargs):
        d = {"vecx": vec_x_func, "vecy": vec_y_func}
        super(PerpendicularVectorFieldDirectionOperator, self).__init__(d, C.none, **kwargs)

    def operation(self, maps):
        mvx = maps["vecx"]
        mvy = maps["vecy"]
        map = arctan2(mvy, mvx)  # in [-numpy.pi, numpy.pi]
        return map


class MassWeightedDensityOperator(FractionOperator):
    r"""
    Mass-weighted density operator

    Parameters
    ----------
    rho_field: ``string``
        name of the density field.
    rho_unit: :class:`~pymses.utils.constants.unit.Unit`
        density field base unit.

    Examples
    --------
    >>> op = MassWeightedDensityOperator("rho", ro.info["unit_density"])

    .. math::
        \left< \rho \right> _{\rho} = \frac{\int\limits_{V} \rho^{2} \mathrm{d}V}{\int\limits_{V} \rho \mathrm{d}V}
    """

    def __init__(self, rho_field, rho_unit):
        if not isinstance(rho_field, str):
            raise AttributeError("'rho_field' attribute must be a valid density field name (string).")
        if not isinstance(rho_unit, C.Unit):
            raise AttributeError("'rho_unit' attribute must be a valid pymses.utils.constants.unit.Unit instance.")

        def num_func(dset):
            if isinstance(dset, IsotropicExtPointDataset):
                return dset[rho_field] ** 2 * dset.get_sizes() ** 3
            else:
                return dset[rho_field] ** 2

        def denom_func(dset):
            if isinstance(dset, IsotropicExtPointDataset):
                return dset[rho_field] * dset.get_sizes() ** 3
            else:
                return dset[rho_field]

        super(MassWeightedDensityOperator, self).__init__(num_func, denom_func, rho_unit)


class MassWeightedGenericScalarOperator(FractionOperator):
    r"""
    Generic mass-weighted scalar field operator

    Parameters
    ----------
    generic_field: ``string``
        name of a generic scalar field (X).
    rho_field: ``string``
        name of the density field.
    generic_field_unit: :class:`~pymses.utils.constants.unit.Unit`
        generic scalar field base unit.

    Examples
    --------
    >>> op = MassWeightedGenericScalarOperator("P", "rho", ro.info["unit_pressure"])

    .. math::
        \left< P \right> _{\rho} = \frac{\int\limits_{V} P \times \rho \mathrm{d}V}{\int\limits_{V} \rho \mathrm{d}V}
    """

    def __init__(self, generic_field, rho_field, generic_field_unit):
        if not isinstance(generic_field, str) or not isinstance(rho_field, str):
            raise AttributeError("'generic_field' and 'rho_field' attributes must be valid field names (string).")
        if not isinstance(generic_field_unit, C.Unit):
            raise AttributeError(
                "'generic_field_unit' attribute must be a valid pymses.utils.constants.unit.Unit instance.")

        def num_func(dset):
            if isinstance(dset, IsotropicExtPointDataset):
                return dset[rho_field] * dset[generic_field] * dset.get_sizes() ** 3
            else:
                return dset[rho_field] * dset[generic_field]

        def denom_func(dset):
            if isinstance(dset, IsotropicExtPointDataset):
                return dset[rho_field] * dset.get_sizes() ** 3
            else:
                return dset[rho_field]

        super(MassWeightedGenericScalarOperator, self).__init__(num_func, denom_func, generic_field_unit)


class MaxLevelOperator(AbstractOperator):
    r"""
    Max. AMR level of refinement operator class
    """

    def __init__(self):
        super(MaxLevelOperator, self).__init__({"levelmax": lambda dset: 0}, C.none, is_max_alos=True, use_cell_dx=True)

    def operation(self, int_dict):
        map = list(int_dict.values())[0]
        return map


class MinimumTemperatureOperator(AbstractOperator):
    """
    Parameters
    ----------
    pressure_field: ``string``
        name of the pressure field.
    rho_field: ``string``
        name of the density field.
    temperature_unit: :class:`~pymses.utils.constants.unit.Unit`
        temperature field unit.

    Examples
    --------
    >>> tempmin_op = MinimumTemperatureOperator("P", "rho", ro.info["unit_temperature"])

    .. math::
        T_{\textrm{min}} = \textrm{min}_{\textrm{line-of-sight}} \left( \frac{P}{\rho} \right)
    """

    def __init__(self, pressure_field, density_field, temperature_unit):
        if not isinstance(pressure_field, str) or not isinstance(density_field, str):
            raise AttributeError("'pressure_field' and 'density_field' attributes must be valid field names (string).")
        if not isinstance(temperature_unit, C.Unit):
            raise AttributeError(
                "'temperature_unit' attribute must be a valid pymses.utils.constants.unit.Unit instance.")

        def invT_func(dset):
            P = dset[pressure_field]
            rho = dset[density_field]
            T = zeros_like(P)
            msk = (rho > 0.0) * (P > 0.0)
            T[msk] = rho[msk] / P[msk]
            T[~msk] = 0.0
            return T

        super(MinimumTemperatureOperator, self).__init__({"invTemp": invT_func}, temperature_unit, is_max_alos=True)

    def operation(self, int_dict):
        return 1.0 / int_dict["invTemp"]


__all__ = ["ScalarOperator", "FractionOperator", "MassWeightedDensityOperator", "MassWeightedGenericScalarOperator",
           "MaxLevelOperator", "MultiFieldOperator", "MinimumTemperatureOperator",
           "PerpendicularVectorFieldDirectionOperator"]
