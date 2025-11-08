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
PyMSES package

"""
from .hpc_conf import HighPerformanceComputingConfiguration
from .etc.rcsetup import rcConfiguration
from .sources.ramses.output import RamsesOutput


# HPC config.
hpc_config = HighPerformanceComputingConfiguration()

# ---------------------------------------------- PyMSES configuration -------------------------------------------------#
rcConfig = rcConfiguration()
del rcConfiguration
# ---------------------------------------------------------------------------------------------------------------------#


# Those following informations are automatically filled by ../update_version_number.py script
# from ../setup.py version number and ../.hg/tags.cache informations
__version__  = '4.1.5'
__revision__ = '$Revision: 1258$'
__date__     = '$Date: 2017-05-05 $'


__all__ = ["RamsesOutput"]
