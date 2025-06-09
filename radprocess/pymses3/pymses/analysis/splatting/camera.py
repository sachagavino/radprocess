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
:mod:`pymses.analysis.splatting.camera` --- Custom splatting camera module
--------------------------------------------------------------------------

"""
import numpy as N
from pymses.utils.regions import Box
from pymses.filters import RegionFilter
from ..camera import Camera
from .convolution_kernels import ConvolKernel


class ExtendedCamera(Camera):
    def __init__(self, camera, ext_box):
        # Camera view area
        filter_box = camera.get_map_box()
        map_range = N.array([[filter_box.min_coords[0], filter_box.max_coords[0]],
                             [filter_box.min_coords[1], filter_box.max_coords[1]],
                             [filter_box.min_coords[2], filter_box.max_coords[2]]])
        far_c = -filter_box.min_coords[2]
        dist = filter_box.max_coords[2]
        rsize = N.diff(map_range, axis=1).transpose()[0]
        # Camera map size + pixel size
        nx_map, ny_map = camera.get_map_size()
        dxy = rsize[:2] / N.array([nx_map, ny_map])

        if isinstance(ext_box, ConvolKernel):
            ext = 3. / 2. * ext_box.max_size
            extension = N.zeros((2, 3))
            ext_axes = ext_box.get_convolved_axes()
            extension[0, ext_axes] = -ext
            extension[1, ext_axes] = ext
        else:  # ext_box is an array
            extension = ext_box

        # Depths limit
        far_c = far_c - extension[0, 2]
        dist = dist + extension[1, 2]

        # Map extension computation : extended map_box
        nplus = (N.ceil(N.abs(extension[:, :2]) / dxy[N.newaxis, :])).astype('i')

        # Unextended map limit coordinates
        ix_min, ix_max = nplus[0, 0] + N.array([0, nx_map])
        iy_min, iy_max = nplus[0, 1] + N.array([0, ny_map])

        # Extended max map size
        mms = N.max([nx_map + N.sum(nplus[:, 0]), ny_map + N.sum(nplus[:, 1])])

        nplus[0, :] = -nplus[0, :]
        new_rsize = rsize[:2] + dxy * N.diff(nplus, axis=0)[0]

        super(ExtendedCamera, self).__init__(camera.center, camera.los_axis, camera.up_vector,
                                             region_size=new_rsize, size_unit=camera.size_unit, distance=dist,
                                             far_cut_depth=far_c, map_max_size=mms, log_sensitive=camera.log_sensitive)

        self.ix_min = ix_min
        self.ix_max = ix_max
        self.iy_min = iy_min
        self.iy_max = iy_max

    def get_window_mask(self):
        """
        Returns the unextended map mask of the camera
        """
        return (self.ix_min, self.ix_max, self.iy_min, self.iy_max)


class CameraFilter(RegionFilter):
    def __init__(self, source, camera, ext_size=0.0, ext_3D=True, use_camera_lvlmax=True):
        """
        Filter build to fit the camera bounding box

        Parameters
        ----------
        source   : pymses data source

        camera   : pymses camera
                        position (along w axis) of the fixed point in the right eye rotation
        ext_size	: float (default 0.0)
            extension to the camera box (in box unit between 0 and 1), used only if ext_3D==True
        ext_3D :	boolean (default True)
            if true, use an ExtendedCamera to extend the filtered region
        use_camera_lvlmax : ``boolean`` (default True)
            Limit the transformation of the AMR grid to particles to AMR cells under the camera octree levelmax
            (so that visible cells are only the ones that have bigger size than the camera pixel size).
        """
        self.cam = camera
        self.cam_box = camera.get_map_box()
        self.ext_3D = ext_3D
        if ext_3D:
            # Extended camera
            r = N.ones((2, 3)) * ext_size
            r[0, :] = -r[0, :]
            ecam = ExtendedCamera(camera, r)
        else:
            ecam = camera

        # Init Filter with extended camera bounding box region
        super(CameraFilter, self).__init__(ecam.get_bounding_box(), source)

        # Max. read level setup
        if use_camera_lvlmax:
            lreq = min(camera.get_required_resolution(), source.get_read_levelmax())
            self.set_read_levelmax(lreq)

    def filtered_dset(self, dset):
        xform = self.cam.viewing_angle_transformation()
        rot_points = xform.transform_points(dset.points)
        if self.ext_3D:
            mask = N.zeros(dset.npoints, dtype=bool)
            sizes = dset.get_sizes()
            unique_sizes = N.unique(sizes)
            if len(unique_sizes) > 4:
                # limit extension size as too big (neighbouring?) cells
                # are normally less interesting as the user zoom on more refined levels
                max_ext_size = unique_sizes[4]
            elif len(unique_sizes) > 0:
                max_ext_size = unique_sizes[-1]
            for size in unique_sizes:
                if size > max_ext_size:
                    ext_size = max_ext_size
                else:
                    ext_size = size
                minb = self.cam_box.min_coords - ext_size
                maxb = self.cam_box.max_coords + ext_size
                b = Box([minb, maxb])
                mask = mask + (b.contains(rot_points) * (size == sizes))
        else:
            mask = self.cam_box.contains(rot_points)

        return dset.filtered_by_mask(mask)


__all__ = ["CameraFilter", "ExtendedCamera"]

