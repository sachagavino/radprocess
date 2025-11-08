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
:mod:`pymses.analysis.cube3d.cube` --- 3D datacube module
---------------------------------------------------------

"""
import os
import numpy as N

from pymses.utils.filename import FileUtil
from ..camera import Camera
from pymses.utils import constants as C, HDF5Serializable
from pymses.utils.constants.unit import Unit


class Datacube3D(HDF5Serializable):
    """
    3D datacube object

    Parameters
    ----------
    camera: :class:`~pymses.analysis.camera.Camera`
        view camera
    cube_size: ``float``
        cube size. Must be strictly positive.
    nside: ``int``
        number of voxels along each cube edge (must be > 2).
    unit: :class:`~pymses.utils.constants.unit.Unit`
        datacube values unit.
    name: ``string``
        datacube name. Default None.
    """

    def __init__(self, camera, cube_size, nside, unit, name=None):
        if not isinstance(camera, Camera):
            raise AttributeError("camera is not a valid Camera instance.")
        self._cam = camera

        if not isinstance(cube_size, float) or cube_size <= 0.0:
            raise AttributeError("'cube_size' attribute must be a positive float value.")
        self._cube_size = cube_size

        if nside < 2:
            raise AttributeError("Not enough voxels along each cube edge : must be grater or equal than 2.")
        self._nside = nside

        if not isinstance(unit, C.Unit):
            raise AttributeError("'unit' attribute is not a valid Unit instance.")
        else:
            self._unit = unit

        if name is not None and not isinstance(name, str):
            raise AttributeError("'name' attribute, if set, must be a string")
        self.name = name

        self._datacube = N.zeros((nside, nside, nside), dtype='d', order='C')

    def __eq__(self, other):
        if not N.allclose(self._datacube, other.data, rtol=1.0e-6):
            return False

        if not self._cam != other.camera:
            return False

        if self.name != other.name:
            return False
        return True

    @property
    def data(self):
        return self._datacube

    @property
    def camera(self):
        return self._cam

    @property
    def unit(self):
        return self._unit

    @property
    def nside(self):
        return self._nside

    @property
    def size(self):
        return self._cube_size

    @property
    def edge_coordinates(self):
        return N.linspace(self._cube_size/2., self._cube_size/2., self._nside + 1)

    def save_fits(self, fits_fname, axis_unit=None):
        """
        Save the Datacube3D into a FITS file

        fits_fname: ``string``
            FITS file name
        axis_unit: : class:`~pymses.utils.constants.unit.Unit`
            output scale unit. Default None (use the camera size unit).
        """
        try:
            from astropy import wcs
            from astropy.io import fits
        except ImportError:
            raise ImportError("astropy package is not available. It is mandatory to save FITS files.")

        cam = Camera(self._cam.center, line_of_sight_axis=self._cam.los_axis, up_vector=self._cam.up_vector,
                     region_size=[self._cube_size, self._cube_size], distance=self._cube_size/2.,
                     far_cut_depth=self._cube_size/2., map_max_size=self._nside, size_unit=self._cam.size_unit)

        hdr = cam.get_fits_header(nz=self._nside, axis_unit=axis_unit)
        hdr["BUNIT"] = self._unit.name

        # Save into FITS file.
        cube = N.transpose(self._datacube, axes=(2, 1, 0))
        hdu = fits.PrimaryHDU(cube, header=hdr)
        hdulist = fits.HDUList([hdu])

        # Save PPVCube, overwriting the FITS file if it already exists.
        path = FileUtil.new_filepath(fits_fname, append_extension=FileUtil.FITS_FILE)
        hdulist.writeto(path, clobber=True)
        print("FITS File '%s' saved." % path)

    def _h5_serialize(self, h5group, **kwargs):
        """
        Serialize the Datacube3D  object into a HDF5 file.

        Parameters
        ----------
        h5group: ``h5py.Group``
            Main group to serialize the Datacube3D into.
        float32: ``bool``
            save datacube values with dtype=float32 instead of float64 ? default False.
        """
        float32 = kwargs.get("float32", False)

        # Save camera object
        cam_group = h5group.create_group("camera")
        self._cam.save_HDF5(cam_group)

        if float32:
            # use appropriate type to save disk space :
            cub = self._datacube.astype("float32")
        else:
            cub = self._datacube

        cub_group = h5group.create_group("datacube")

        # Save cube name, if any
        if self.name is not None:
            cub_group.attrs["name"] = str(self.name)

        # Save data
        cub_group.create_dataset("data", data=cub)

        # Save meta information
        cub_group.attrs["cube_size"] = self.size
        cub_group.attrs["nside"] = self.nside

        # Save datacube value unit
        unit_group = cub_group.create_group("unit")
        self._unit.save_HDF5(unit_group)

    @classmethod
    def _h5_deserialize(cls, h5group, version):
        """
        Read a Datacube3D object from a HDF5 file.

        Parameters
        ----------
        h5group: ``h5py.Group``
            Main group to deserialize the Datacube3D from.
        version: ``int``
            Version of the Datacube3D class to deserialize.

        Returns
        -------
        cube: :class:`~pymses.analysis.cube3d.cube.Datacube3D`
            new Datacube3D instance.
        """
        if version == Datacube3D._pymses_h5_version:  # Current Datacube3D version
            cam = Camera.from_HDF5(h5group["camera"])
            cube_group = h5group['datacube']

            if "name" in cube_group.attrs:
                name = cube_group.attrs['name']
            else:
                name = None

            cube_size = cube_group.attrs['cube_size']
            nside = cube_group.attrs['nv']

            unit = Unit.from_HDF5(cube_group["unit"])

            c = cls(cam, cube_size, nside, unit, name=name)
            c.data[...] = cube_group['data'][...]
        else:
            raise ValueError("Unknown Datacube3D version (%d)" % int(version))

        return c


__all__ = ["Datacube3D"]
