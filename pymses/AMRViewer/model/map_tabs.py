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
from pymses.analysis import Camera
from numpy import array, newaxis
from .ramses import map_engine, MapEngine
from .region_finder import RegionFinderModel
from .utils import get_appropriate_unit


class MapTabsModel(object):
    def __init__(self):
        self.unitDict = {}
        self.mapUnit = {}
        self.hdf5_map_file = None
        self.refresh_densityRT_tab = False
        self.refresh_levelmax_tab = False
        self.map_levelmax = None

    def reset(self):
        self.tab_map_arrays = {}
        self.glass_x = None
        self.glass_y = None
        self.rule_dist = None
        self.hdf5_map_file = None
        self.refresh_densityRT_tab = False
        self.refresh_levelmax_tab = False
        self.map_levelmax = None

    def UpdateTabMap(self, nxy, tab_name):
        cam = self.get_los_camera(nxy)
        cam.size_unit = self.ro.info["unit_length"]
        if self.rf_lower_resolution is not None:
            cam.map_max_size = self.tab_lower_resolution
        if self.source is not None and (tab_name == "levelmax" or tab_name == "densityRT"):
            # Special case : compute two images tab with the same map processor
            if "densityRT" not in self.map_engine_dict:
                self.map_engine_dict["densityRT"] = map_engine(self, "densityRT", self.ro, self.source)
            me = self.map_engine_dict["densityRT"]
            if "levelmax" not in self.map_engine_dict:
                self.map_engine_dict["levelmax"] = map_engine(self, "levelmax", self.ro, self.source)
            me_levelmax = self.map_engine_dict["levelmax"]
            self.map_name_list["densityRT"], self.map_name_list["levelmax"], self.map_levelmax = \
                me.MakeBothMaps({"los_zoom": cam}, me_levelmax)
            self.tab_map_arrays["densityRT"] = me.GetMapArrays()["los_zoom"]
            self.unitDict["densityRT"] = me.getUnitDict()
            self.mapUnit["densityRT"] = me.getMapUnit()
            self.image_computed["densityRT"] = True

            self.tab_map_arrays["levelmax"] = me_levelmax.GetMapArrays()["los_zoom"]
            self.unitDict["levelmax"] = me_levelmax.getUnitDict()
            self.mapUnit["levelmax"] = me_levelmax.getMapUnit()
            self.image_computed["levelmax"] = True
            if tab_name == "levelmax":
                self.refresh_densityRT_tab = True
                print("DensityRT tab computed too !")
            else:
                self.refresh_levelmax_tab = True
                print("Levelmax tab computed too !")
        else:
            # Normal case : compute image for one tab only
            if tab_name not in self.map_engine_dict:
                self.map_engine_dict[tab_name] = map_engine(self, tab_name, self.ro, self.source)
            me = self.map_engine_dict[tab_name]
            if tab_name == "stars_particles":
                cache_dset = self.stars_particles_files_cache_dset
            elif tab_name == "dm_particles":
                cache_dset = self.dm_particles_files_cache_dset
            else:
                cache_dset = self.amr_files_cache_dset
            self.map_name_list[tab_name] = me.MakeMaps({"los_zoom": cam}, self.cmap, False, self.FFTkernelSizeFactor)
            if not self.rememberSomeData:
                cache_dset.clear()
            m = me.GetMapArrays()
            self.tab_map_arrays[tab_name] = m["los_zoom"]
            self.unitDict[tab_name] = me.getUnitDict()
            self.mapUnit[tab_name] = me.getMapUnit()
            self.image_computed[tab_name] = True

    def GetTabImage(self, tab_name, verbose=True):
        if tab_name not in self.map_engine_dict:
            self.map_engine_dict[tab_name] = map_engine(self, tab_name, self.ro, self.source)
        me = self.map_engine_dict[tab_name]
        if self.adaptive_gaussian_blur and self.map_levelmax is not None:
            return me.MakeImage(self.map_name_list[tab_name], ["los_zoom"], self.cmap, \
                                self.map_levelmax, ramses_output=self.ro, verbose=verbose)["los_zoom"]
        return me.MakeImage(self.map_name_list[tab_name], ["los_zoom"], self.cmap, \
                            self.adaptive_gaussian_blur, ramses_output=self.ro, verbose=verbose)["los_zoom"]

    def SetGlassCenter(self, position):
        """
        Sets the position of the magnifying glass center
        """
        self.glass_x, self.glass_y = position

    def GetMap(self, tab_name):
        return self.tab_map_arrays.get(tab_name)

    def GetMapLengthUnit(self):
        return self.length_unit

    def GetGlassValue(self, tab_name):
        map = self.GetMap(tab_name)
        if ((self.glass_x is None) + (self.glass_y is None)):
            return (None, None)
        else:
            if self.glass_x >= map.shape[0] or self.glass_y >= map.shape[1] or self.glass_x < 0 or self.glass_y < 0:
                return (None, None)
            val = map[self.glass_x, self.glass_y]
            v = val * self.mapUnit.get(tab_name)
            return get_appropriate_unit(v, self.unitDict[tab_name])

    def SetRuleDist(self, dist):
        if self.rule_dist == dist:
            return False
        else:
            self.rule_dist = dist
            return True

    def GetRuleValue(self):
        if self.rule_dist is None:
            return (0.0, "")
        else:
            size = self.rule_dist * self.region_size * self.length_unit
            if self.zoom_size is not None:
                size = size * self.zoom_size
            return get_appropriate_unit(size, RegionFinderModel.lunit_dict)

    def LoadmapFromHDF5(self, h5fname):
        import tables as T
        tab_name = "hdf5"
        print("Openning file : ", h5fname)
        self.map_name_list[tab_name] = [h5fname]
        h5f = T.openFile(h5fname, 'r')
        cam = Camera.from_HDF5(h5f)
        if cam.log_sensitive:
            map = 10 ** h5f.getNode("/map/map").read()
        else:
            map = h5f.getNode("/map/map").read()
        # FLIP_TOP_BOTTOM :
        l = map.shape[1] - 1
        for j in range(l + 1):
            map[:, l - j] = map[:, j]
        self.tab_map_arrays[tab_name] = map
        try:
            # TODO: try to get proper unit:
            self.mapUnit[tab_name] = h5f.getNode("/map/unit/dimensions").read() * h5f.getNode("/map/unit/val").read()
            self.unitDict[tab_name] = {self.mapUnit[tab_name]: self.mapUnit[tab_name]}
            print(h5f.getNode("/map/unit/dimensions").read())
            print(h5f.getNode("/map/unit/val").read())
        except:
            self.mapUnit[tab_name] = C.none
            self.unitDict[tab_name] = {"(?)": C.none}
        self.image_computed[tab_name] = True
        h5f.close()

        #self.LoadLosCameraFromHDF5(h5fname)

    def get_los_camera(self, nxy):
        uaxis, vaxis = self.GetUVAxes()
        waxis = self.GetLosAxis()
        cc = self.region_center
        rs = self.region_size
        cam = Camera(center=cc, line_of_sight_axis=waxis, up_vector=vaxis, region_size=[rs, rs],
                     distance=(rs / 2.), far_cut_depth=(rs / 2.), map_max_size=nxy)
        return cam

