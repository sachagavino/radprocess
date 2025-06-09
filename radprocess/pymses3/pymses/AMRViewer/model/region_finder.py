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

from numpy import array, newaxis
from .ramses import *
from ..view.dialogs import *
from .utils import get_appropriate_unit


class RegionFinderModel(object):
    lunit_dict = {'pc': C.pc,
                  'kpc': C.kpc,
                  'Mpc': C.Mpc,
                  'Gpc': C.Gpc,
                  'ly': C.ly,
                  'au': C.au}

    def __init__(self):
        self.length_unit = None
        self.map_name_list_history = []

    def reset(self):
        print("Reset Region finder Model")
        self.cross_center = array([0.5, 0.5, 0.5])
        self.region_center = array([0.5, 0.5, 0.5])
        self.cam_dict = {'los': None, 'u': None, 'v': None}

        self.cursor_u = None
        self.cursor_v = None
        self.cursor_w = None

        self.SetRegionSize(1.0)
        self.SetZoomSize(None)

    def region_finder_previous_cam(self):
        if self.map_name_list["regionFinder"] == self.map_name_list_history[-1]:
            self.map_name_list_history.pop()
        self.map_name_list["regionFinder"] = self.map_name_list_history[-1]
        self.LoadLosCameraFromFile(self.map_name_list["regionFinder"][0])

    def UpdateCursorUVWPosition(self, u, v, w):
        """
        Update the (u, v, w) position of the cursor.
        """
        cs = self.region_size
        if u is not None:
            self.cursor_u = u * cs
        else:
            self.cursor_u = None
        if v is not None:
            self.cursor_v = v * cs
        else:
            self.cursor_v = None
        if w is not None:
            self.cursor_w = w * cs
        else:
            self.cursor_w = None

        self.update_cross_center_coords()

    def UpdateRegionFinderMaps(self, nxys):
        self.update_cameras_dict(nxys)
        if self.regionFinderMapType == "stars_particles":
            cache_dset = self.stars_particles_files_cache_dset
        elif self.regionFinderMapType == "dm_particles":
            cache_dset = self.dm_particles_files_cache_dset
        else:
            cache_dset = self.amr_files_cache_dset
            if self.source_camera is not None and not self.source_camera.contains_camera(self.cam_dict["los"]):
                dlg = ConfirmBuildSourceDialog()
                dlg.ShowModal()
                if (dlg.answer == 1):
                    # Yes button
                    cam = self.cam_dict["los"].copy()
                    # extend loading by 2 AMR level:
                    cam.map_max_size = cam.map_max_size * 4
                    self.BuildSource(cam)
                elif (dlg.answer < 0):
                    # Cancel button
                    return False
                dlg.Destroy()
        if self.regionFinderMapType not in self.map_engine_dict:
            self.map_engine_dict[self.regionFinderMapType] = map_engine(self,
                                                                        self.regionFinderMapType, self.ro, self.source)
        me = self.map_engine_dict[self.regionFinderMapType]
        if self.rf_lower_resolution is not None:
            self.cam_dict["los"].map_max_size = self.rf_lower_resolution
            self.cam_dict["u"].map_max_size = self.rf_lower_resolution
            self.cam_dict["v"].map_max_size = self.rf_lower_resolution
        self.map_name_list["regionFinder"] = me.MakeMaps(self.cam_dict, self.cmap, self.multiprocessing,
                                                         self.FFTkernelSizeFactor)
        self.map_name_list_history.append(self.map_name_list["regionFinder"])
        if not self.rememberSomeData:
            cache_dset.clear()
        self.image_computed["regionFinder"] = True
        return True

    def UpdateImage(self, verbose=True):
        if self.regionFinderMapType not in self.map_engine_dict:
            self.map_engine_dict[self.regionFinderMapType] = map_engine(self,
                                                                        self.regionFinderMapType, self.ro, self.source)
        me = self.map_engine_dict[self.regionFinderMapType]
        self.axis_images = me.MakeImage(self.map_name_list["regionFinder"],
                                        list(self.cam_dict.keys()), self.cmap, self.adaptive_gaussian_blur,
                                        ramses_output=self.ro, verbose=verbose)

    def GetAxisImages(self):
        return self.axis_images

    def GetCrossCenterUVW(self):
        u = self.cursor_u
        if u is None:
            u = 0.0
        v = self.cursor_v
        if v is None:
            v = 0.0
        w = self.cursor_w
        if w is None:
            w = 0.0
        uvw_center = array([u, v, w])
        return uvw_center

    def update_cross_center_coords(self):
        uvw_center = self.GetCrossCenterUVW()
        self.cross_center = self.cam_dict["los"].deproject_points(uvw_center[newaxis, :])[0]

    def update_cameras_dict(self, nxys):
        uaxis, vaxis = self.GetUVAxes()
        waxis = self.GetLosAxis()

        self.region_center[:] = self.cross_center[:]
        cc = self.region_center

        if (self.zoom_size is not None):
            self.region_size = self.region_size * self.zoom_size
        rs = self.region_size
        self.cam_dict['los'] = Camera(center=cc, line_of_sight_axis=waxis, up_vector=vaxis, region_size=[rs, rs],
                                      size_unit=self.length_unit, distance=(rs / 2.), far_cut_depth=(rs / 2.),
                                      map_max_size=nxys["los"])
        self.cam_dict['u'] = Camera(center=cc, line_of_sight_axis=uaxis, up_vector=vaxis, region_size=[rs, rs],
                                    size_unit=self.length_unit, distance=(rs / 2.), far_cut_depth=(rs / 2.),
                                    map_max_size=nxys["u"])
        self.cam_dict['v'] = Camera(center=cc, line_of_sight_axis=-vaxis, up_vector=waxis, region_size=[rs, rs],
                                    size_unit=self.length_unit, distance=(rs / 2.), far_cut_depth=(rs / 2.),
                                    map_max_size=nxys["v"])

    def SaveLosCameraToFile(self, file_path):
        cam = self.get_los_camera(256)
        fname = basename(file_path)
        ext = fname[-3:]
        if (ext == ".h5" or ext == ".hdf5"):
            cam.save_HDF5(file_path)
        else:
            cam.save_csv(file_path)
        return True

    def LoadLosCameraFromFile(self, file_path):
        fname = basename(file_path)
        ext = fname[-3:]
        if (ext == ".h5" or ext == ".hdf5"):
            self.LoadLosCamera(Camera.from_HDF5(file_path))
        else:
            self.LoadLosCamera(Camera.from_csv(file_path))

    def LoadLosCamera(self, cam):
        self.SetRegionSize(cam.region_size[0])
        uvect, vvect, losvect = cam.get_camera_axis()
        self.SetUVAxis(u=uvect, v=vvect)
        self.RefreshPoints()
        self.region_center = cam.center
        self.cam_dict["los"] = cam
        self.update_cross_center_coords()

    def SetRegionSize(self, rsize):
        self.region_size = rsize

    def SetRegionSizeUnit(self, unit):
        self.length_unit = unit

    def SetZoomSize(self, s):
        self.zoom_size = s

    def GetCrossCenter(self):
        return self.cross_center

    def GetRegionSize(self):
        if self.ro is None:
            return (0.0, "")
        else:
            size = self.region_size * self.length_unit
            if self.zoom_size is not None:
                size = size * self.zoom_size
            return get_appropriate_unit(size, RegionFinderModel.lunit_dict)
