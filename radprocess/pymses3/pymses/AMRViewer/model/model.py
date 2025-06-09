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

from .ramses import RamsesModel
from .axes3d import Axes3DModel
from .region_finder import RegionFinderModel
from .map_tabs import MapTabsModel


class Model(RamsesModel, Axes3DModel, RegionFinderModel, MapTabsModel):
    def __init__(self):
        RamsesModel.__init__(self)
        Axes3DModel.__init__(self)
        RegionFinderModel.__init__(self)
        MapTabsModel.__init__(self)
        self.map_engine_dict = {}
        self.image_computed = {}
        self.amr_files_cache_dset = {}
        self.stars_particles_files_cache_dset = {}
        self.dm_particles_files_cache_dset = {}
        self.map_name_list = {}

    def reset(self):
        print("Model reset")
        Axes3DModel.reset(self)
        RegionFinderModel.reset(self)
        MapTabsModel.reset(self)
        RamsesModel.reset(self)
        self.map_engine_dict.clear()
        self.image_computed.clear()
        self.amr_files_cache_dset.clear()
        self.stars_particles_files_cache_dset.clear()
        self.dm_particles_files_cache_dset.clear()
        self.map_name_list.clear()

    def freeTheMemory(self):
        self.amr_files_cache_dset.clear()
        self.stars_particles_files_cache_dset.clear()
        self.dm_particles_files_cache_dset.clear()
        RamsesModel.reset(self)
