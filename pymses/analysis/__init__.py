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
:mod:`pymses.analysis` --- Analysis and post-processing package
###############################################################

"""
from .camera import Camera, SphericalCamera
from .datamap import DataMap
from .operator import FractionOperator, ScalarOperator, MultiFieldOperator, MaxLevelOperator,\
    MassWeightedGenericScalarOperator, MassWeightedDensityOperator, MinimumTemperatureOperator, \
    PerpendicularVectorFieldDirectionOperator


from .profile_binners import *
from .avg_point import *

__all__ = ['Camera', 'SphericalCamera',
           'DataMap', 'MultiFieldOperator', 'ScalarOperator', 'FractionOperator', 'MaxLevelOperator',
           'MassWeightedGenericScalarOperator', 'MassWeightedDensityOperator', 'MinimumTemperatureOperator',
           'PerpendicularVectorFieldDirectionOperator']
