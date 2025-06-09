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
:mod:`pymses.analysis.camera` --- Camera module
-----------------------------------------------

"""
from collections import OrderedDict
import numpy as N

from pymses import __version__ as pymses_ver
from pymses.core.transformations import *
from pymses.utils import HDF5Serializable, constants as C
from pymses.utils.regions import Box
from pymses.utils._point_utils import meshgrid
from pymses.utils._ray_cube_utils import clip_unit_size_cube as _clip_cube
from .transfer_functions import ColorLinesTransferFunction


class Camera(HDF5Serializable):
    r"""
    Camera class for 2D projected maps computing. Take a look at documentation figures to get a clearer definition.

    Parameters
    ----------
    center: ``list`` of 3 ``float`` values
        region of interest center coordinates (default value is [0.5, 0.5, 0.5], the simulation domain center).
    line_of_sight_axis :
        axis of the line of sight (z axis is the default_value) : not zero-normed [ux, uy, uz] array or simulation
        domain specific axis key "x", "y" or "z".
    up_vector          :
        direction of the y axis of the camera (up). If None, the up vector is set to the z axis (or y axis if the
        line-of-sight is set to the z axis). If given a not zero-normed [ux, uy, uz] array is expected (or a simulation
        domain specific axis key "x", "y" or "z").
    region_size:
        projected size of the region of interest (default [1.0, 1.0])
    size_unit          : :class:`~pymses.utils.constants.unit.Unit`
        length unit (default : 1 astronomical unit).
    distance           :
        distance of the camera from the center of interest (along the line-of-sight axis, default 0.5).
    far_cut_depth      :
        distance of the background (far) cut plane from the center of interest (default 0.5). The region of interest
        is within the camera position and the far cut plane.
    map_max_size       :
        maximal resolution of the camera (default 1024 pixels)
    log_sensitive: ``bool``
        whether the camera pixels are log sensitive or not (default True).

    Examples
    --------
        >>> cam  = Camera(center=[0.5, 0.5, 0.5], line_of_sight_axis='z', up_vector='y', region_size=[1., 1.],
        ... size_unit=ramses_output.info["unit_length"], distance=0.5, far_cut_depth=0.5, map_max_size=1024)

    """
    special_axes = OrderedDict([("x", N.array([1., 0., 0.])),
                                ("y", N.array([0., 1., 0.])),
                                ("z", N.array([0., 0., 1.]))])

    def __init__(self, center=None, line_of_sight_axis="z", up_vector=None, region_size=None,
                 size_unit=C.pc, distance=0.5, far_cut_depth=0.5, map_max_size=512, log_sensitive=True):
        super(Camera, self).__init__()

        # Region of interest center coordinates
        if center is None:
            self.center = N.array([0.5, 0.5, 0.5])
        else:
            self.center = N.asarray(center, dtype='d')

        # Axis of the line of sight
        self.los_axis = self._get_input_cam_axis(line_of_sight_axis)

        # Size unit
        self.size_unit = size_unit

        # Up vector of the camera
        if up_vector is not None:
            self.up_vector = self._get_input_cam_axis(up_vector)
            if N.allclose(N.cross(self.up_vector, self.los_axis), 0.0, rtol=1.0e-2):
                raise AttributeError("Up vector must be different from line-of-sight axis")
        else:
            # Default up vector of the camera : z axis. if the line of sight is the z (or -z) axis,
            # the up vector is set to y-axis
            if N.isclose(N.abs(N.sum(self.los_axis * Camera.special_axes["z"])), 1.0):
                self.up_vector = Camera.special_axes["y"]
            else:
                self.up_vector = Camera.special_axes["z"]

        # Map extent
        if region_size is not None:
            self.region_size = N.asarray(region_size, 'd')
        else:
            self.region_size = N.array([1., 1.])

        # Region of interest depth
        self.distance = distance
        self.far_cut_depth = far_cut_depth

        # Map maximal size (pixels)
        self.map_max_size = map_max_size

        # Logarithmic sensitivity of the Camera captors
        self.log_sensitive = log_sensitive

        # Set the color transfer function to None
        self.color_tf = None

    def _get_input_cam_axis(self, axis):
        """
        Handles input axis attribute 'x', '-x', 'y', '-y', 'z', '-z' or any 3D axis.

        Parameters
        ----------
        axis: axis attribute

        Returns
        -------
        v: 3D axis
        """
        if isinstance(axis, str):
            if len(axis) == 1 and axis in self.special_axes:
                return self.special_axes[axis]
            elif len(axis) == 2 and axis[0] == '-' and axis[1] in self.special_axes:
                return -self.special_axes[axis[1]]
            else:
                raise AttributeError("Unknown axis '%s" % axis)
        else:
            return axis / N.linalg.norm(axis, 2)

    def same_axes(self, other_cam):
        if not isinstance(other_cam, Camera):
            raise AttributeError("'other_cam' must be a Camera instance. Got %s." % type(other_cam))
        axes = self.get_camera_axis()
        oaxes = other_cam.get_camera_axis()
        return N.allclose(axes, oaxes, rtol=1.0e-6)

    def __eq__(self, other):
        """
        Check whether a camera is equal to an other Camera

        Parameters
        ----------
        other: :class:`~pymses.analysis.camera.Camera`
             Other Camera instance to compare to.

        Returns
        -------
        eq: ``bool``
            True if center, line-of-sight axes, map boxes and resolution of both Camera are
            identical, otherwise False.
        """
        if not isinstance(other, Camera):
            return False

        if not N.allclose(self.center, other.center, rtol=1.0e-6):
            return False

        if not self.same_axes(other):
            return False

        bmin, bmax = self.get_map_box().get_bounding_box()
        obmin, obmax = other.get_map_box().get_bounding_box()
        if not N.allclose(bmin, obmin, rtol=1.0e-6) or not N.allclose(bmax, obmax, rtol=1.0e-6):
            return False

        # Same size unit ?
        if self.size_unit != other.size_unit:
            return False

        if self.map_max_size != other.map_max_size:
            return False

        return True

    def get_axis_names(self):
        """
        Get the list of the camera axis names

        Returns
        -------
        l : ``list`` of ``string``
            list of the camera axis names
        """
        axes = self.get_camera_axis()
        default = ["u", "v", "w"]
        names = []
        for iaxis in range(axes.shape[0]):
            axis = axes[iaxis, :]
            n = None
            for saxis_name, saxis in Camera.special_axes.items():
                if N.allclose(saxis, axis, rtol=1.0e-6):
                    n = saxis_name
                    break
                elif N.allclose(-saxis, axis, rtol=1.0e-6):
                    n = "-%s" % saxis_name
                    break
            if n is not None:
                names.append(n)
            else:
                names.append(default[iaxis])

        return names

    def get_3D_cameras(self, z_fixed_point=0.0, perspective_angle_deg=2.0, depth_shift_percent=0.0):
        """
        Get the 3D left/right eye cameras for stereoscopic views.

        Parameters
        ----------

        perspective_angle_deg: ``float``
            perspective angle between left and right eye camera (in degrees, default 2.0).
        z_fixed_point: ``float``
            position (along line-of-sight axis) of the fixed point in the left/right eye perspective. Default 0 (center
            of the camera).
        depth_shift_percent: ``float``
            global depth shift (in % of the camera width), positive values shifts the camera box content forward, and
            negative values shifts it backward. Default is 0 (no depth shift).

        Returns
        -------

        (left_camera, right_camera): the left/right eye Camera object tuple for 3D image processing
        """
        if perspective_angle_deg < 0.0:
            raise AttributeError("perspective angle must pe a positive or null float value !")

        if perspective_angle_deg == 0.0 and depth_shift_percent == 0.0:
            raise AttributeError("either perspective angle or depth shift must be set to a non null value.")

        if depth_shift_percent <= -100.0 or depth_shift_percent >= 100.0:
            raise AttributeError("depth shift must be expressed in percent (in the range ]-100.0, 100[).")
        # Left/right eye camera
        lcam = self.copy()
        rcam = self.copy()

        if perspective_angle_deg > 0.0:
            if z_fixed_point != 0.0:
                lcam.center = self.center + self.los_axis * z_fixed_point
                rcam.center = lcam.center

            # Rotate left/right Camera around vaxis
            # and keep up_vector to vaxis (Important! prevent spurious effects when los ~ up_vector)
            vaxis = self.get_camera_axis()[1]
            lcam = lcam.get_rotated_camera(ang_deg=-perspective_angle_deg / 2.0, use_up_vector=False)
            lcam.up_vector = vaxis[:]
            rcam = rcam.get_rotated_camera(ang_deg=perspective_angle_deg / 2.0, use_up_vector=False)
            rcam.up_vector = vaxis[:]


            if z_fixed_point != 0.0:
                lcam.center = lcam.center - lcam.los_axis * z_fixed_point
                rcam.center = rcam.center - rcam.los_axis * z_fixed_point

        if depth_shift_percent != 0.0:
            # Depth-shifting the left/right cameras
            horiz_shift = self.region_size[0] * depth_shift_percent / 200.0
            luaxis = lcam.get_camera_axis()[0, :]
            ruaxis = lcam.get_camera_axis()[0, :]
            lcam.center -= horiz_shift * luaxis
            rcam.center += horiz_shift * ruaxis

        cams = (lcam, rcam)
        return cams

    def get_required_resolution(self):
        """
        Returns
        -------

        lev : ``int``
            the level of refinement up to which one needs to read the data to compute the projection
            of the region of interest with the specified resolution.

        """
        lev = int(N.ceil(N.log2(self.map_max_size / max(self.region_size))))
        return lev

    def get_map_size(self):
        """
        Returns
        -------
        (nx, ny) : (``int``, ``int``) ``tuple``
            the size (nx,ny) of the image taken by the camera (pixels)

        """
        aspect_ratio = float(self.region_size[0]) / float(self.region_size[1])
        if aspect_ratio > 1.:
            nx_map = self.map_max_size
            ny_map = int(N.round(self.map_max_size / aspect_ratio))
        else:
            nx_map = int(N.round(self.map_max_size * aspect_ratio))
            ny_map = self.map_max_size

        return nx_map, ny_map

    def contains_camera(self, cam):
        """
        Parameters
        ----------
        An other camera object

        Returns
        ----------
        Boolean : True if data needed for this camera view include all data
            needed for the camera view given in argument.

        """
        return (self.get_required_resolution() >= cam.get_required_resolution()) and \
               self.get_bounding_box().contains(cam.get_bounding_box())

    def get_map_box(self):
        """
        Returns the (0.,0.,0.) centered cubic bounding box of the area covered by the camera
        """
        dx, dy = self.region_size
        bound_min = N.array([-dx / 2., -dy / 2., -self.far_cut_depth], dtype='d')
        bound_max = N.array([dx / 2., dy / 2., self.distance], dtype='d')
        box = Box([bound_min, bound_max])
        return box

    def get_rotated_camera(self, ang_deg=1.0, use_up_vector=False):
        """
        TODO
        """
        # if ang_deg <= 0.0:
        #     raise AttributeError("'ang_deg' parameter must be > 0.0")

        if use_up_vector:
            j = self.up_vector
        else:
            uaxis, vaxis, los = self.get_camera_axis()
            j = vaxis

        R = rot3d_axvector(j, ang_deg * N.pi / 180.)
        cam = self.copy()
        cam.los_axis = R.transform_points(N.array([self.los_axis]))[0]
        return cam

    def get_pixel_surface(self):
        """ Returns the surface of any pixel of the camera
        """
        dx, dy = self.region_size
        nx, ny = self.get_map_size()
        S = 1.0 * dx / nx * dy / ny
        return S

    def get_pixels_coordinates_centers(self):
        """
        Returns the center values of the camera pixels x/y coordinates. The pixel coordinates of the center of the
        camera is the origin (0,0).
        """
        xedges, yedges = self.get_pixels_coordinates_edges()
        xcen = (xedges[1:] + xedges[:-1]) / 2.
        ycen = (yedges[1:] + yedges[:-1]) / 2.
        centers = (xcen, ycen)
        return centers

    def get_pixels_coordinates_edges(self):
        """
        Returns the edges values of the camera pixels x/y coordinates
        The pixel coordinates of the center of the camera is the origin (0,0)
        """
        dx, dy = self.region_size
        nx, ny = self.get_map_size()
        xedges = N.linspace(-float(dx) / 2., float(dx) / 2., num=nx + 1)
        yedges = N.linspace(-float(dy) / 2., float(dy) / 2., num=ny + 1)

        edges = (xedges, yedges)
        return edges

    def get_uvaxes_edges_labels(self, axis_unit=None, force_zero_centered_axis_values=False):
        """
        Get the labels and pixel edge coordinates of the camera u and v axes.

        Parameters
        ----------
        axis_unit: :class:`pymses.utils.constants.unit.Unit` or ``None``
            user-defined axis unit (default: ``None``). If left to ``None``, the camera size unit is used.
        force_zero_centered_axis_values: ``bool``
            Used when any map axis matches one of the x/y/z axes. If False, displays the x/y/z axis coordinate value.
             Otherwise displays the 0-centered u/v axis coordinate value. default False.

        Returns
        -------
        labedges: ``tuple``of ``tuple``
            ((u_axisname, u_axisunit, u_label_latex, uedges), (v_axisname, v_axisunit, v_label_latex, vedges)) tuple of
            the label and pixel edge coordinates of the camera u and v axes.
        """
        if axis_unit is not None and not isinstance(axis_unit, C.Unit):
            raise AttributeError("'axis_unit' must be a valid Unit instance.")

        uedges, vedges = self.get_pixels_coordinates_edges()
        uaxis, vaxis, zaxis = self.get_camera_axis()

        # U-axis (pointing upward)
        u_axisname = 'u'
        iaxis = 0
        for axname, axis in Camera.special_axes.items():  # Try to match the u/v axes with the x/y/z axes
            prod = N.sum(uaxis * axis)
            if N.isclose(N.abs(prod), 1.0, rtol=1.0e-6):
                u_axisname = axname
                if not force_zero_centered_axis_values:
                    uedges += self.center[iaxis]
                if prod < 0.0:
                    uedges = uedges[::-1]
                break
            iaxis += 1

        # V-axis (pointing rightward)
        v_axisname = 'v'
        iaxis = 0
        for axname, axis in Camera.special_axes.items():  # Try to match the u/v axes with the x/y/z axes
            prod = N.sum(vaxis * axis)
            if N.isclose(N.abs(prod), 1.0, rtol=1.0e-6):
                v_axisname = axname
                if not force_zero_centered_axis_values:
                    vedges += self.center[iaxis]
                if prod < 0.0:
                    vedges = vedges[::-1]
                break
            iaxis += 1

        # Handle size unit => update edge coordinates + axis labels
        if axis_unit is not None:
            axisunit = axis_unit.name
            u_label_latex = "%s (%s)" % (u_axisname, axis_unit.latex)
            v_label_latex = "%s (%s)" % (v_axisname, axis_unit.latex)
            size_fact = self.size_unit.express(axis_unit)
            uedges *= size_fact
            vedges *= size_fact
        else:
            axisunit = self.size_unit.name
            u_label_latex = "%s (%s)" % (u_axisname, self.size_unit.latex)
            v_label_latex = "%s (%s)" % (v_axisname, self.size_unit.latex)

        t = ((u_axisname, axisunit, u_label_latex, uedges), (v_axisname, axisunit, v_label_latex, vedges))
        return t

    def get_fits_header(self, nz=1, axis_unit=None):
        """
        Get the FITS file header information from this Camera instance attributes.

        Parameters
        ----------
        nz: ``int``
            number of Z-axis slices (1: 2D, nz>1: 3D). Default 1 (2D).
        axis_unit:
            TODO

        Returns
        -------
        header: ``dict``
            FITS file header dictionary
        """
        try:
            from astropy import wcs
        except ImportError:
            raise ImportError("astropy module is not available")

        nx, ny = self.get_map_size()

        # Get u/v axes labels and pixel edge coordinates
        uinfo, vinfo = self.get_uvaxes_edges_labels(axis_unit=axis_unit)
        _dummy1, axisunit, u_label_latex, uedges = uinfo
        _dummy2, axisunit, v_label_latex, vedges = vinfo
        uname, vname, wname = self.get_axis_names()

        dx = uedges[1]-uedges[0]
        dy = vedges[1]-vedges[0]

        try:
            if nz == 1:  # 2D map
                w = wcs.WCS(naxis=2)
                w.wcs.crpix = [nx / 2, ny / 2]  # Map center
                w.wcs.cdelt = [dx, dy]  # Coordinate increment
                w.wcs.crval = [0, 0]
                w.wcs.ctype = [uname, vname]
                w.wcs.cunit = 2 * ['%s' % axisunit]

                hdr = w.to_header()

                hdr['BITPIX'] = -32
                hdr['NAXIS'] = 2
                hdr['NAXIS1'] = nx
                hdr['NAXIS2'] = ny

            elif nz < 1:
                raise AttributeError("nz must be a strictly positive integer.")
            else:  # (nz > 1) 3D datacube
                dz = (self.distance + self.far_cut_depth) / nz
                nzc = int(self.distance / dz)
                if axis_unit is not None:
                    dz *= self.size_unit.express(axis_unit)

                w = wcs.WCS(naxis=3)
                w.wcs.crpix = [nx / 2, ny / 2, nzc]  # Map center
                w.wcs.cdelt = N.array([dx, dy, dz])  # Coordinate increment
                w.wcs.crval = [0, 0, 0]
                w.wcs.ctype = [uname, vname, wname]
                w.wcs.cunit = 3 * ['%s' % axisunit]

                hdr = w.to_header()
                hdr['BITPIX'] = -32
                hdr['NAXIS'] = 3
                hdr['NAXIS1'] = nx
                hdr['NAXIS2'] = ny
                hdr['NAXIS3'] = nz

            hdr['ORIGIN'] = 'RAMSES'
            hdr['CREATOR'] = 'PyMSES %s' % pymses_ver
            hdr['XCENTER'] = (self.center[0], "Camera center x-axis coordinate")
            hdr['YCENTER'] = (self.center[1], "Camera center y-axis coordinate")
            hdr['ZCENTER'] = (self.center[2], "Camera center z-axis coordinate")

            u, v, los = self.get_camera_axis()
            hdr['ABS_X'] = (u[0], "X-axis coordinate of the map horizontal axis")
            hdr['ABS_Y'] = (u[1], "Y-axis coordinate of the map horizontal axis")
            hdr['ABS_Z'] = (u[2], "Z-axis coordinate of the map horizontal axis")
            hdr['ORD_X'] = (v[0], "X-axis coordinate of the map vertical axis")
            hdr['ORD_Y'] = (v[1], "Y-axis coordinate of the map vertical axis")
            hdr['ORD_Z'] = (v[2], "Z-axis coordinate of the map vertical axis")
            hdr['LOS_X'] = (los[0], "Line-of-sight x-axis coordinate")
            hdr['LOS_Y'] = (los[1], "Line-of-sight y-axis coordinate")
            hdr['LOS_Z'] = (los[2], "Line-of-sight z-axis coordinate")
            return hdr
        except Exception as ex:
            print("Cannot generate FITS file header from this Camera instance.")
            raise ex

    def viewing_angle_rotation(self):
        """Returns the rotation corresponding to the viewing angle of the camera
        """
        cam_axis = self.get_camera_axis()
        rot = LinearTransformation(cam_axis)
        return rot

    def viewing_angle_transformation(self):
        """Returns the transformation corresponding to the viewing angle
        of the camera
        """
        rot = self.viewing_angle_rotation()
        tr = translation(-self.center)
        Tr = rot * tr
        return Tr

    def get_camera_axis(self):
        """Returns the camera u, v and z axis coordinates
        """
        z_axis = self.los_axis
        u_axis = N.cross(self.up_vector, self.los_axis)
        u_axis = u_axis / N.linalg.norm(u_axis)
        v_axis = N.cross(z_axis, u_axis)
        cam_axis = N.array([u_axis, v_axis, z_axis])
        return cam_axis

    def get_region_size_level(self):
        """ Returns the level of the AMR grid for which
        the cell size ~ the region size
        """
        lev = int(N.round(N.log2(1. / max(self.region_size))))
        return lev

    def get_slice_points(self, z=0.0):
        """
        Returns the (x, y, z) coordinates of the points contained in a slice plane
        perpendicular to the line-of-sight axis at a given position z.

        z --- slice plane position along line-of-sight (default 0.0 => center of the region)
        """
        if not isinstance(z, float):
            raise AttributeError("'z' attribute must be a float value.")

        uc, vc = self.get_pixels_coordinates_centers()
        points = meshgrid([uc, vc, N.array([z], dtype='d')])
        return self.deproject_points(points)

    def get_datacube_info(self, cube_size, resolution=256):
        """
        Returns the (x, y, z) coordinate array of the points contained in a 3D datacube and the required max. AMR level
        to read.

        Parameters
        ----------
        cube_size: ``float``
            3D datacube size (must be positive)
        resolution: `` int``
            datacube number of voxels along each dimension. Default 256.

        Returns
        -------
        points: ``numpy.ndarray``
            (npoints, 3) shaped point coordinate array.
        rlev: ``int``
            Required AMR level to read.
        """
        if not isinstance(cube_size, float) or cube_size <= 0.0:
            raise AttributeError("'cube_size' attribute must be a positive float value.")

        edges = N.linspace(-float(cube_size) / 2., float(cube_size) / 2., num=resolution + 1)
        cellcen = (edges[1:] + edges[:-1]) / 2.
        uvw_points = meshgrid([cellcen, cellcen, cellcen])

        # Required AMR level to read
        rlev = int(N.ceil(N.log2(resolution / cube_size)))
        points = self.deproject_points(uvw_points)
        return points, rlev

    def get_bounding_box(self):
        """Returns the bounding box of the region of interest in the simulation
        domain corresponding of the area covered by the camera
        """
        b = self.get_map_box()
        box_bounds = [N.array([b.min_coords[0], b.max_coords[0]], dtype='d'),
                      N.array([b.min_coords[1], b.max_coords[1]], dtype='d'),
                      N.array([b.min_coords[2], b.max_coords[2]], dtype='d')]
        corner_list = meshgrid(box_bounds)
        xform_corners = self.deproject_points(corner_list)
        xform_min_bounds = N.min(xform_corners, axis=0)
        xform_max_bounds = N.max(xform_corners, axis=0)
        min_bounds = N.max([xform_min_bounds, N.zeros(3)], axis=0)
        max_bounds = N.min([xform_max_bounds, N.ones(3)], axis=0)
        box = Box([min_bounds, max_bounds])
        return box

    def deproject_points(self, uvw_points, origins=None):
        """
        Return xyz_coords deprojected coordinates of a set of points from given [u,v,w] coordinates

        Parameters
        ----------
        - (u=0,v=0, w=0) is the center of the camera.
        - v is the coordinate along the vaxis
        - w is the depth coordinate of the points along the	line-of-sight of the camera.
        if origins is True, perform a vectorial transformation of the vectors described by uvw_points
        anchored at positions 'origins'
        """
        xform = self.viewing_angle_transformation().inverse()
        if origins is None:
            coords_xyz = xform.transform_points(uvw_points)
        else:
            coords_xyz = xform.transform_vectors(uvw_points, origins)
        return coords_xyz

    def project_points(self, points):
        """
        Return a (coords_uv, depth) tuple where 'coord_uv' is the projected coordinates of
        a set of points on the camera plane. (u=0,v=0) is the center of the camera plane.
        'depth' is the depth coordinate of the points along the line-of-sight of the camera.

        Parameters
        ----------
        points				: ``numpy array of floats``
                array of points(x,y,z) coordinates to project
        """
        xform = self.viewing_angle_transformation()
        xform_pts = xform.transform_points(points)
        coords_uv = xform_pts[:, :2]
        depth = xform_pts[:, 2]
        return (coords_uv, depth)

    def get_map_mask(self, float32=True):
        """
        Returns the mask map of the camera. Each pixel has an alpha parameter between :
        * 1, if the ray of the pixel intersects the simulation domain (unit size cube) along its entire length,
        * and 0, if it does not intersects the simulation domain.

        Parameters
        ----------
        float32		``Boolean`` (default True)
            use float32 numpy dtype array instead of float64 to save memory
        """
        # Camera map size
        nx_map, ny_map = self.get_map_size()

        # Camera bounding box
        b = self.get_map_box()
        len_max = b.max_coords[2] - b.min_coords[2]
        box_bounds = [N.array([b.min_coords[0], b.max_coords[0]], dtype='d'),
                      N.array([b.min_coords[1], b.max_coords[1]], dtype='d'),
                      N.array([b.min_coords[2], b.max_coords[2]], dtype='d')]
        corner_list = meshgrid(list(box_bounds))
        xform_corners = self.deproject_points(corner_list)

        # The region of interest is completely inside the simulation domain
        if (N.min(xform_corners) >= 0.).all() and (N.max(xform_corners) <= 1.0).all():
            mask = N.ones((nx_map, ny_map))
            if float32:
                mask = mask.astype("float32")
            return mask
        # Get camera rays
        _dummy_vects, _dummy_orig, ray_length = self.get_rays(clip_unit_size_cube=True)

        mask = ray_length
        if len_max > 0.:
            mask /= len_max

        mask = mask.reshape(nx_map, ny_map)
        if float32:
            mask = mask.astype("float32")
        return mask

    def get_rays(self, clip_unit_size_cube=True, custom_origins=None):
        """
        Returns ray_vectors, ray_origins and ray_lengths arrays for ray tracing ray definition
        """
        # Camera view area
        centered_map_box = self.get_map_box()

        # Min./max. camera depth
        xmin, ymin, zmin = centered_map_box.min_coords[:]
        xmax, ymax, zmax = centered_map_box.max_coords[:]

        if custom_origins is not None:
            if custom_origins.shape[1] != 3:
                raise AttributeError("'custom_origins' must be a 3D point array.")

            # Total number of rays
            n_rays = custom_origins.shape[0]

            # custom origin points projected in the camera system of coordinates
            ray_orig = N.zeros((n_rays, 3))
            ray_orig[:, :2], ray_orig[:, 2] = self.project_points(custom_origins)

            # Camera background plane clipping
            msk = ray_orig[:, 2] < zmin
            ray_orig[msk, 2] = zmin
            msk = (ray_orig[:, 0] >= xmin) * (ray_orig[:, 0] <= xmax)
            msk *= (ray_orig[:, 1] >= ymin) * (ray_orig[:, 1] <= ymax)
            msk *= ray_orig[:, 2] <= zmax
            ray_lengths = N.zeros(n_rays)
            ray_lengths[msk] = zmax - ray_orig[msk, 2]
        else:
            # Camera map size
            nx_map, ny_map = self.get_map_size()

            # Total number of rays
            n_rays = nx_map*ny_map

            # Pixels edges coordinates
            xcen, ycen = self.get_pixels_coordinates_centers()

            # Ray origins => grid on the camera background cut plane
            ray_orig = N.zeros((nx_map, ny_map, 3))
            ray_orig[:, :, 0] = xcen[:, N.newaxis]
            ray_orig[:, :, 1] = ycen[N.newaxis, :]
            ray_orig[:, :, 2] = zmin
            ray_orig = ray_orig.reshape(n_rays, 3)

            ray_lengths = N.ones(n_rays) * (zmax - zmin)

        # Origin points (in AMR grid coordinates : x, y, z)
        ray_origins = self.deproject_points(ray_orig)

        # Axis of the rays (in AMR grid coordinates : x, y, z)
        ray_vectors = N.zeros((n_rays, 3))
        ray_vectors[:, 2] = 1.0
        ray_vectors = self.deproject_points(ray_vectors, ray_origins)

        # Simulation domain ([0; 1]^3 cube) clipping, if required
        if clip_unit_size_cube:
            ray_ends = ray_origins + ray_lengths[:, N.newaxis] * ray_vectors
            if N.min(ray_origins) < 0.0 or N.max(ray_origins) > 1.0 or N.min(ray_ends) < 0.0 or N.max(ray_ends) > 1.0:
                _clip_cube(ray_origins, ray_vectors, ray_lengths)

        return ray_vectors, ray_origins, ray_lengths

    def set_color_transfer_function(self, tf):
        assert isinstance(tf, ColorLinesTransferFunction)
        self.color_tf = tf

    def _h5_serialize(self, h5group, **kwargs):
        """
        Serialize the camera parameters into a HDF5 file.

        Parameters
        ----------
        h5group: ``h5py.Group``
            HDF5 group object where to store the camera parameters.
        """
        b = self.get_map_box()
        bound_uv_depth = N.array([b.min_coords, b.max_coords])
        h5group.create_dataset("center", data=self.center)
        h5group.create_dataset("axis_coordinates", data=self.get_camera_axis())
        h5group.create_dataset("bound_uv_depth", data=bound_uv_depth)
        size_unit_group = h5group.create_group("size_unit")
        self.size_unit.save_HDF5(size_unit_group)
        h5group.attrs['map_max_size'] = self.map_max_size
        h5group.attrs['log_sensitive'] = self.log_sensitive

    @classmethod
    def _h5_deserialize(cls, h5group, version):
        """
        Read a camera from a HDF5 file.

        Parameters
        ----------
        h5group: ``h5py.Group``
            Main group to deserialize the Camera from.
        version: ``int``
            Version of the Camera class to deserialize.

        Returns
        -------
        cam: :class:`~pymses.analysis.camera.Camera`
            new Camera instance
        """
        # Check Camera object version written in HDF5 file
        if version == cls._pymses_h5_version:  # Current camera version (h5py syntax)
            center = h5group['center'][:]
            uaxis, vaxis, los = h5group['axis_coordinates'][...]
            bounds = h5group['bound_uv_depth'][...]
            fcd = -bounds[0, 2]
            dist = bounds[1, 2]
            region_size = N.diff(bounds, axis=0)[0, :-1]
            su = C.Unit.from_HDF5(h5group['size_unit'])
            mms = h5group.attrs['map_max_size']
            is_log_sens = h5group.attrs['log_sensitive']
        else:
            raise AttributeError("Unknown Camera version (%d)" % int(version))

        cam = cls(center=center, line_of_sight_axis=los, up_vector=vaxis, region_size=region_size, size_unit=su,
                  distance=dist, far_cut_depth=fcd, map_max_size=mms, log_sensitive=is_log_sens)
        return cam

    @classmethod
    def _load_legacy_cam(cls, cam_group):
        """
        Read  a camera from a legacy (PyMSES v<=4.0.0) HDF5 group.

        Parameters
        ----------
        cam_group: :class:`h5py.Group`
            Main group to deserialize the Camera from

        Returns
        -------
        cam: :class:`~pymses.analysis.camera.Camera`
            camera instance
        """
        # Old camera version (PyTables syntax)
        center = cam_group.center.read()
        uaxis, vaxis, los = cam_group.axis_coordinates.read()
        bounds = cam_group.bound_uv_depth.read()
        fcd = -bounds[0, 2]
        dist = bounds[1, 2]
        region_size = N.diff(bounds, axis=0)[0, :-1]
        su = C.pc  # Default size unit : parsec
        mms = cam_group.map_max_size.read()
        is_log_sens = cam_group.log_sensitive.read()
        return cls(center=center, line_of_sight_axis=los, up_vector=vaxis, region_size=region_size, size_unit=su,
                   distance=dist, far_cut_depth=fcd, map_max_size=mms, log_sensitive=is_log_sens)

    @classmethod
    def load_legacy_HDF5_camera(cls, h5fname):
        """
        Read  a camera from a legacy (PyMSES v<=4.0.0) HDF5 file.

        Parameters
        ----------
        h5fname: ``string``
            HDF5 file name

        Returns
        -------
        cam: :class:`~pymses.analysis.camera.Camera`
            camera instance
        """
        import h5py
        try:
            import tables as T
        except ImportError:
            raise ImportError("PyTables module is required to handle legacy HDF5 files in PyMSES !")

        ######################################################################################
        #      Old Pymses (v<=4.0.0) HDF5 camera/map files written with PyTables check       #
        ######################################################################################
        f0 = h5py.File(h5fname, 'r')
        if 'camera' not in f0 or'CAMERA_VERSION' in f0['camera'].attrs:
            f0.close()
            raise AttributeError("'%s' file is not a valid (legacy) Camera File" % h5fname)
        f0.close()
        ######################################################################################
        f = T.File(h5fname, 'r')
        g = f.get_node("/camera")
        try:
            cam = cls._load_legacy_cam(g)
        except Exception as exc:
            raise IOError("HDF5 I/O error : %s" % exc)
        finally:
            f.close()
        return cam

    def copy(self):
        """
        Returns a copy of this camera.
        """
        cen = self.center.copy()
        los = self.los_axis.copy()
        vaxis = self.up_vector.copy()
        rs = self.region_size.copy()
        return Camera(center=cen, line_of_sight_axis=los, up_vector=vaxis, region_size=rs, size_unit=self.size_unit,
                      distance=self.distance, far_cut_depth=self.far_cut_depth, map_max_size=self.map_max_size,
                      log_sensitive=self.log_sensitive)

    def __unicode__(self):
        """
        Format Camera object as a string
        """
        return "[Camera]: center = %s, line of sight = %s,\n" \
               "          up_vector = %s, region_size= %s,\n" \
               "          distance = %g, far_cut_depth = %g, map_max_size = %d, log_sensitive = %s" %\
               (self.center, self.los_axis, self.up_vector, self.region_size, self.distance, self.far_cut_depth,
                self.map_max_size, self.log_sensitive)


class SphericalCamera(Camera):
    """
    Spherical projection camera subclass used for virtual reality (VR) image production
    """
    _default_image_width_to_height_ratio = 2.0
    _default_stereo_delta_angle_degrees = 0.0
    def __init__(self, *args, **kwargs):
        delta_deg = kwargs.pop('stereo_delta_deg', SphericalCamera._default_stereo_delta_angle_degrees)
        self._img_wh_ratio = kwargs.pop('image_wh_ratio', SphericalCamera._default_image_width_to_height_ratio)
        self._delta_rad = delta_deg/180.0 * N.pi

        # Front plane clipping
        # zenith_clfactor = 8.  # 50.0
        # horizon_clfactor = 5.  # 10.0
        horz_clip_angle_deg = kwargs.pop('horizon_clip_angle_deg', 6.0)  # 5.0 * IPD_2
        self._alpha_clip_horz_rad = horz_clip_angle_deg / 180.0 * N.pi
        zen_clip_angle_deg = kwargs.pop('zenith_clip_angle_deg', 3.5)  # 8.0 * IPD_2
        self._alpha_clip_zen_rad = zen_clip_angle_deg / 180.0 * N.pi

        super(SphericalCamera, self).__init__(*args, **kwargs)

        b = self.get_map_box()
        box_bounds = [N.array([b.min_coords[0], b.max_coords[0]], dtype='d'),
                      N.array([b.min_coords[1], b.max_coords[1]], dtype='d'),
                      N.array([b.min_coords[2], b.max_coords[2]], dtype='d')]
        corner_list = meshgrid(box_bounds)
        self._radius = N.max(N.sqrt(N.sum(corner_list ** 2, axis=1)))

        corner_list = meshgrid([N.array([0.0, 1.0], dtype='d'),
                                N.array([0.0, 1.0], dtype='d'),
                                N.array([0.0, 1.0], dtype='d')])
        max_radius = N.max(N.sqrt(N.sum((corner_list-self.center[N.newaxis, :])** 2, axis=1)))
        if self._radius > max_radius:
            self._radius = max_radius

    @property
    def radius(self):
        return self._radius

    @property
    def image_wh_ratio(self):
        return self._img_wh_ratio

    @image_wh_ratio.setter
    def image_wh_ratio(self, r):
        self._img_wh_ratio = r

    @property
    def stereo_delta_deg(self):
        return self._delta_rad / N.pi * 180.0

    @stereo_delta_deg.setter
    def stereo_delta_deg(self, d):
        self._delta_rad = d / 180.0 * N.pi

    @property
    def horizon_clip_angle_deg(self):
        return self._alpha_clip_horz_rad / N.pi * 180.0

    @horizon_clip_angle_deg.setter
    def horizon_clip_angle_deg(self, d):
        self._alpha_clip_horz_rad = d / 180.0 * N.pi

    @property
    def zenith_clip_angle_deg(self):
        return self._alpha_clip_zen_rad / N.pi * 180.0

    @zenith_clip_angle_deg.setter
    def zenith_clip_angle_deg(self, d):
        self._alpha_clip_zen_rad = d / 180.0 * N.pi

    def get_map_mask(self, float32=True):
        """
        Returns the mask map of the camera. Each pixel has an alpha parameter between :
        * 1, if the ray of the pixel intersects the simulation domain (unit size cube) along its entire length,
        * and 0, if it does not intersects the simulation domain.

        Parameters
        ----------
        float32		``Boolean`` (default True)
            use float32 numpy dtype array instead of float64 to save memory
        """
        # Camera map size
        nx_map, ny_map = self.get_map_size()

        # Get camera rays
        _dummy_vects, _dummy_orig, ray_length, _dS, _front_clip = self.get_rays(clip_unit_size_cube=True)

        mask = N.ones_like(ray_length)
        mask[N.isclose(ray_length, 0.0)] = 0.0

        mask = mask.reshape(nx_map, ny_map)
        if float32:
            mask = mask.astype("float32")
        return mask

    def get_bounding_box(self):
        """Returns the bounding box of the region of interest in the simulation
        domain corresponding of the area covered by the camera
        """
        min_bounds = N.max([self.center - self._radius, N.zeros(3)], axis=0)
        max_bounds = N.min([self.center + self._radius, N.ones(3)], axis=0)
        return Box([min_bounds, max_bounds])

    def get_map_size(self):
        """
        Returns
        -------
        (nx, ny) : (``int``, ``int``) ``tuple``
            the size (nx,ny) of the image taken by the camera (pixels)

        """
        nlongitude = self.map_max_size

        # Map paramters
        nlat = int(N.round(nlongitude / self._img_wh_ratio))

        return nlongitude, nlat

    def get_rays(self, clip_unit_size_cube=True, custom_origins=None):
        """
        Returns ray_vectors, ray_origins and ray_lengths + pixel solid angle arrays for ray tracing ray definition
        Also returns the minimum
        """
        nlongitude, nlat = self.get_map_size()

        # Total number of rays
        n_rays = nlongitude * nlat

        # longitude/latitude coordinate values
        phi_e = N.linspace(-N.pi/2., N.pi/2., nlat+1)
        lam_e = N.linspace(-N.pi, N.pi, nlongitude+1)
        phi = (phi_e[1:] + phi_e[:-1]) / 2.
        lam = (lam_e[1:] + lam_e[:-1]) / 2.

        # Ray vectors
        ray_vects = N.zeros((nlongitude, nlat, 3))
        ray_vects[:, :, 0] = -N.cos(phi[N.newaxis, :])*N.sin(lam[:, N.newaxis])
        ray_vects[:, :, 1] = -N.sin(phi[N.newaxis, :])
        ray_vects[:, :, 2] = N.cos(phi[N.newaxis, :])*N.cos(lam[:, N.newaxis])
        ray_vectors = ray_vects.reshape(n_rays, 3)
        ray_vectors = self.deproject_points(ray_vectors, self.center[N.newaxis, :])

        # Ray origins/ends/lengths
        ray_origins = self.center[N.newaxis, :] - self._radius * ray_vectors

        # delta = -0.05  # < 0 (left), > 0.0 (right)
        if self._delta_rad != 0.0:
            IPD_2 = self._radius * N.abs(N.sin(self._delta_rad))  # Half of the interpupillary distance
            side_sign = 1
            if self._delta_rad < 0.0:
                side_sign = -1
            print("IPD = %g" % (IPD_2 * 2.0))
            ray_ends = N.zeros((nlongitude, nlat, 3))
            ray_ends[:, :, 0] = N.cos(lam[:, N.newaxis]) * IPD_2 * side_sign
            ray_ends[:, :, 1] = 0.0
            ray_ends[:, :, 2] = N.sin(lam[:, N.newaxis]) * IPD_2 * side_sign
            ray_ends = ray_ends.reshape(n_rays, 3)
            ray_ends = self.deproject_points(ray_ends)
            ray_vectors = ray_ends - ray_origins
            ray_vectors /= N.linalg.norm(ray_vectors, ord=2, axis=-1)[:, N.newaxis]
            horizon_clip_dist = IPD_2 / N.tan(self._alpha_clip_horz_rad)
            zenith_clip_dist = IPD_2 / N.tan(self._alpha_clip_zen_rad)
            front_clip = N.ones_like(lam)[:, N.newaxis] * (zenith_clip_dist - (zenith_clip_dist - horizon_clip_dist) *
                                                           N.cos(phi[N.newaxis, :]))
            front_clip = front_clip.reshape(n_rays)
            ray_lengths = N.ones(n_rays) * self._radius * (1.0 + N.sin(self._delta_rad)**2)
        else:
            ray_lengths = N.ones(n_rays) * self._radius
            front_clip = N.zeros_like(ray_lengths, dtype='d')
            ray_ends = N.ones_like(ray_origins) * self.center[N.newaxis, :]

        # Simulation domain ([0; 1]^3 cube) clipping, if required
        if clip_unit_size_cube:
            if N.min(ray_origins) < 0.0 or N.max(ray_origins) > 1.0 or N.min(ray_ends) < 0.0 or N.max(ray_ends) > 1.0:
                _clip_cube(ray_origins, ray_vectors, ray_lengths)

        # Pixel solid angle
        dcosphi = N.diff(N.sin(phi_e))
        dlam = N.diff(lam_e)
        dS = dcosphi[N.newaxis, :] * dlam[:, N.newaxis]  # steradian
        dS = dS.reshape(n_rays)

        return ray_vectors, ray_origins, ray_lengths, dS, front_clip

    def get_3D_cameras(self, perspective_angle_deg=2.0):
        """
        Get the 3D left/right eye cameras for stereoscopic views.

        Parameters
        ----------

        perspective_angle_deg: ``float``
            perspective angle between left and right eye camera (in degrees, default 2.0).

        Returns
        -------
        (left_camera, right_camera): the left/right eye Camera object tuple for 3D image processing
        """
        if perspective_angle_deg < 0.0:
            raise AttributeError("perspective angle must pe a positive or null float value !")

        # Left/right eye camera
        lcam = self.copy()
        rcam = self.copy()

        if perspective_angle_deg > 0.0:
            lcam.stereo_delta_deg = -perspective_angle_deg / 2.0
            rcam.stereo_delta_deg = perspective_angle_deg / 2.0

        cams = (lcam, rcam)
        return cams

    def copy(self):
        """
        Returns a copy of this spherical camera.
        """
        cen = self.center.copy()
        los = self.los_axis.copy()
        vaxis = self.up_vector.copy()
        rs = self.region_size.copy()
        return SphericalCamera(center=cen, line_of_sight_axis=los, up_vector=vaxis, region_size=rs,
                               size_unit=self.size_unit, distance=self.distance, far_cut_depth=self.far_cut_depth,
                               map_max_size=self.map_max_size, log_sensitive=self.log_sensitive,
                               stereo_delta_deg=self.stereo_delta_deg, image_wh_ratio=self.image_wh_ratio,
                               horizon_clip_angle_deg=self.horizon_clip_angle_deg,
                               zenith_clip_angle_deg=self.zenith_clip_angle_deg)

    def __eq__(self, other):
        """
        Check whether a spherical camera object is equal to an other SphericalCamera instance

        Parameters
        ----------
        other: :class:`~pymses.analysis.camera.SphericalCamera`
             Other SphericalCamera instance to compare to.

        Returns
        -------
        eq: ``bool``
            True if both SphericalCamera are identical, otherwise False.
        """
        if not isinstance(other, SphericalCamera):
            return False

        if not super(SphericalCamera, self).__eq__(other):
            return False

        return self.stereo_delta_deg == other.stereo_delta_deg

    def _h5_serialize(self, h5group, **kwargs):
        """
        Serialize the spherical camera parameters into a HDF5 file.

        Parameters
        ----------
        h5group: ``h5py.Group``
            HDF5 group object where to store the spherical camera parameters.
        """
        super(SphericalCamera, self)._h5_serialize(h5group, **kwargs)
        h5group.attrs['stero_delta_deg'] = self.stereo_delta_deg
        h5group.attrs['image_wh_ratio'] = self.image_wh_ratio

    @classmethod
    def _h5_deserialize(cls, h5group, version):
        """
        Read a spherical camera from a HDF5 file.

        Parameters
        ----------
        h5group: ``h5py.Group``
            Main group to deserialize the Camera from.
        version: ``int``
            Version of the SphericalCamera class to deserialize.

        Returns
        -------
        cam: :class:`~pymses.analysis.camera.SphericalCamera`
            new SphericalCamera instance
        """
        cam = super(SphericalCamera, cls)._h5_deserialize(h5group, version)

        # Check SphericalCamera object version written in HDF5 file
        if version == cls._pymses_h5_version:  # Current spherical camera version (h5py syntax)
            if 'stero_delta_deg' not in h5group.attrs:
                cam.stereo_delta_deg = SphericalCamera._default_stereo_delta_angle_degrees
            else:
                cam.stereo_delta_deg = h5group.attrs['stero_delta_deg']
            if 'image_wh_ratio' not in h5group.attrs:
                cam.image_wh_ratio = SphericalCamera._default_image_width_to_height_ratio
            else:
                cam.image_wh_ratio = h5group.attrs['image_wh_ratio']

        else:
            raise AttributeError("Unknown Camera version (%d)" % int(version))

        return cam

__all__ = ["Camera", "SphericalCamera"]
