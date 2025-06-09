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
:mod:`pymses.analysis.quantities.base` --- Base quantity module
---------------------------------------------------------------

"""
import numpy as N
from pymses.utils.constants import Unit


class AbstractQuantity(object):
    def __init__(self, name, amr_fields, unit, tex_name=r""):
        self._name = ""
        self.name = name
        self._amr_fields = []
        self.amr_fields = amr_fields
        self._unit = None
        self.unit = unit
        self._tex_name = r""
        self.tex_name = tex_name

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        if not isinstance(new_name, str):
            raise AttributeError("'name' attribute must be a string.")
        self._name = new_name

    @property
    def amr_fields(self):
        return self._amr_fields

    @amr_fields.setter
    def amr_fields(self, new_amr_fields):
        if not isinstance(new_amr_fields, list):
            raise AttributeError("'amr_fields' attribute must be a list.")
        for f in new_amr_fields:
            if not isinstance(f, str):
                raise AttributeError("'amr_fields' attribute must be a list of strings.")

        self._amr_fields = new_amr_fields

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, new_unit):
        if not isinstance(new_unit, Unit):
            raise AttributeError("'unit' attribute must be a valid pymses.constants.Unit instance.")
        self._unit = new_unit

    @property
    def tex_name(self):
        return self._tex_name

    @tex_name.setter
    def tex_name(self, new_tex_name):
        if not isinstance(new_tex_name, str):
            raise AttributeError("'tex_name' attribute must be a string.")
        self._tex_name = new_tex_name

    def func_dset(self, amr_dataset):
        raise NotImplementedError()

    def data_unit(self, info):
        raise NotImplementedError()

    # @classmethod
    # def process(cls, ramses_output, dset=None, unit=None, verbose=False):
    #     if dset == None:
    #         source = ramses_output.amr_source(cls.amrFields,
    #                                grav_compat=GRAV_COMPAT,
    #                                verbose=verbose)
    #         dset = CellsToPoints(source).flatten()
    #     return cls.get_calc(ro, unit)(dset)
    #
    # @classmethod
    # def get_calc(cls, ro, unit=None):
    #     """
    #     Returns a function f evaluating a Pointdataset.
    #     This function f returns : _func(dset)*factor
    #     """
    #     if unit is None:
    #         unit = cls.unit
    #     # varplot = partial(lambda dset,func,factor: func(dset)*factor,
    #     #                      func=cls._func,
    #     #                      factor=cls.factorFunc(ro).express(unit))
    #     varplot = lambda dset: cls._func(dset) * cls.factorFunc(ro).express(unit)
    #     return varplot
    #
    # @classmethod
    # def label(cls):
    #     if cls.unit is C.none:
    #         return r"$%s$" % cls.texName
    #     else:
    #         return r"$%s \ \left( \ %s \ \right)$" %(cls.texName, cls.texUnitName)


class ScalarQuantity(AbstractQuantity):
    @property
    def abs_qty(self):
        class _AbsoluteValueQuantity(self.__class__):
            def func_dset(self, amr_dataset):
                return N.abs(self.func_dset(amr_dataset))

        nm = "{name!s}_abs".format(name=self.name)
        texnm = r"\vert {tex_name!s} \vert".format(tex_name=self.tex_name)
        return _AbsoluteValueQuantity(nm, self.amr_fields, self.unit, texnm)


class VectorQuantity(AbstractQuantity):
    _axes_dict = {'x': 0, 'y': 1, 'z': 2}

    @property
    def tex_name(self):
        return r"\textbf{{{tex_name!s}}}".format(tex_name=self._tex_name)

    @property
    def norm_qty(self):
        class _NormQuantity(AbstractQuantity):
            def func_dset(self, amr_dataset):
                return N.linalg.norm(self.func_dset(amr_dataset), axis=-1)

        nm = "{name!s}_norm".format(name=self.name)
        texnm = r"\Vert {tex_name!s} \Vert".format(tex_name=self.tex_name)
        return _NormQuantity(nm, self.amr_fields, self.unit, texnm)

    def __getitem__(self, item):
        if item in VectorQuantity._axes_dict:
            ax = item
            ind = VectorQuantity._axes_dict[item]
        else:
            for ax, ind in VectorQuantity._axes_dict.items():
                if ind == item:
                    break
            raise KeyError("'{key!s}' item was not found in {class_name!s}".format(key=item,
                                                                                   class_name=self.__class__.__name__))
        class _VectorAxisComponentQty(ScalarQuantity):
            def func_dset(self, amr_dataset):
                dset = self.func_dset(amr_dataset)
                return dset[..., ind]
        nm = "{name!s}_{axis!s}".format(name=self.name, axis=ax)
        texnm = r"{tex_name!s}_{{\textrm{{{axis!s}}}}}".format(tex_name=self._tex_name, axis=ax)
        return _VectorAxisComponentQty(nm, self.amr_fields, self.unit, texnm)

    def dot(self, v):
        class _DotProductQuantity(ScalarQuantity):
            def func_dset(self, amr_dataset):
                return N.dot(self.func_dset(amr_dataset), v) / N.linalg.norm(v)

        nm = "{name!s}_dot_v".format(name=self.name)
        texnm = r"{tex_name!s} \cdot \textbf{{v}}".format(tex_name=self.tex_name)
        return _DotProductQuantity(nm, self.amr_fields, self.unit, texnm)


class MultivaluedQuantity(AbstractQuantity):
    @property
    def sum_qty(self):
        class _SumQuantity(ScalarQuantity):
            def func_dset(self, amr_dataset):
                return N.sum(self.func_dset(amr_dataset), axis=-1)

        nm = "{name!s}_tot".format(name=self.name)
        texnm = r"{tex_name!s}^{{\textrm{{tot}}}}".format(tex_name=self.tex_name)
        return _SumQuantity(nm, self.amr_fields, self.unit, texnm)

    def __getitem__(self, item):
        if not isinstance(item, int):
            raise KeyError("'{key!s}' item was not found in {class_name!s}".format(key=item,
                                                                                   class_name=self.__class__.__name__))
        class _MultivaluedComponentQty(ScalarQuantity):
            def func_dset(self, amr_dataset):
                dset = self.func_dset(amr_dataset)
                return dset[..., item]
        nm = "{name!s}_{item:d}".format(name=self.name, item=item)
        texnm = r"{tex_name!s}_{{{item:d}}}".format(tex_name=self.tex_name, item=item)
        return _MultivaluedComponentQty(nm, self.amr_fields, self.unit, texnm)


__all__ = ['AbstractQuantity', "ScalarQuantity", "VectorQuantity", 'MultivaluedQuantity']
