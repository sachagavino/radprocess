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
:mod:`pymses.analysis.ppv.cube` --- PPV datacube module
-------------------------------------------------------

"""
import os
import numpy as N

from ..camera import Camera
from pymses.utils import constants as C, HDF5Serializable
from pymses.utils.constants.unit import Unit
from pymses.utils.filename import FileUtil


class PPVCube(HDF5Serializable):
    """
    Position-Position-Velocity 3D datacube object

    Parameters
    ----------
    camera: :class:`~pymses.analysis.camera.Camera`
        view camera
    vrange: `list` of `float`
        velocity range [vmin, vmax]. Default [-50.0, 50.0]
    vunit: :class:`~pymses.utils.constants.unit.Unit`
        velocity range unit. Default: km/s
    nv: `int`
        number of velocity bins (must be > 2). Default 100.
    sigma_unit: :class:`~pymses.utils.constants.unit.Unit`
        datacube values surface density unit. Default: mH/cm^2.
    name: `string`
        datacube name. Default None.
    """

    def __init__(self, camera, sigma_unit=None, vrange=None, nv=100, vunit=None, name=None):
        if not isinstance(camera, Camera):
            raise AttributeError("camera is not a valid Camera instance.")
        self._cam = camera

        if sigma_unit is not None:
            self._sigma_unit = sigma_unit
        else:
            self._sigma_unit = C.Unit(name="cm^-2", base_unit=C.mH / C.cm ** 2,
                                      descr="particles per squared centimeter")

        if vrange is None:
            self._vmin = -50.0
            self._vmax = 50.0
        else:
            self._vmin, self._vmax = vrange
            if self._vmin > self._vmax:
                raise AttributeError("vmin < vmax condition is not satisfied in velocity range attribute.")

        if nv < 2:
            raise AttributeError("Not enough velocity bins : must be grater or equal than 2.")
        self._nv = nv

        if vunit is not None:
            self._vunit = vunit
        else:
            self._vunit = C.Unit(name="km/s", base_unit=C.km / C.s, descr="kilometer per second")

        if name is not None and not isinstance(name, str):
            raise AttributeError("'name' attribute, if set, must be a string")
        self.name = name

        self._nx, self._ny = self._cam.get_map_size()
        self._datacube = N.zeros((self._nx, self._ny, self._nv), dtype='d', order='C')

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
    def min_velocity(self):
        return self._vmin

    @property
    def max_velocity(self):
        return self._vmax

    @property
    def velocity_unit(self):
        return self._vunit

    @property
    def unit(self):
        return self._sigma_unit

    @property
    def velocity_nbins(self):
        return self._nv

    @property
    def velocity_bins(self):
        return N.linspace(self._vmin, self._vmax, self._nv + 1)

    def save_fits(self, fits_fname):
        """

        :param fits_fname:
        :return:
        """
        try:
            from astropy import wcs
            from astropy.io import fits
        except ImportError:
            raise ImportError("astropy package is not available. It is mandatory to save FITS files.")

        hdr = self._cam.get_fits_header(nz=self._nv)
        hdr["CRPIX3"] = self._nv / 2  # Midplane velocity index
        hdr["CDELT3"] = (self._vmax - self._vmin) / self._nv  # Velocity coordinate increment
        hdr["CRVAL3"] = (self._vmin + self._vmax) / 2.  # Midplane velocity value]
        hdr["CTYPE3"] = "Velocity"
        hdr["CUNIT3"] = self._vunit.name

        hdr["BUNIT"] = self._sigma_unit.name

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
        Serialize the PPVCube  object into a HDF5 file.

        Parameters
        ----------
        h5group: ``h5py.Group``
            Main group to serialize the PPVCube into.
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

        # Save map data
        nx, ny = self._cam.get_map_size()
        map_size = N.array([nx, ny])
        vel_range = N.array([self._vmin, self._vmax])

        cub_group = h5group.create_group("ppvcube")

        # Save ppv cube name, if any
        if self.name is not None:
            cub_group.attrs["name"] = str(self.name)

        # Save data
        cub_group.create_dataset("data", data=cub)

        # Save meta information
        cub_group.create_dataset("ppsize", data=map_size)
        cub_group.create_dataset("velocity_range", data=vel_range)
        cub_group.attrs["nv"] = self._nv

        # Save units
        Vunit_group = cub_group.create_group("vel_unit")
        self._vunit.save_HDF5(Vunit_group)
        Sunit_group = cub_group.create_group("sigma_unit")
        self._sigma_unit.save_HDF5(Sunit_group)

    @classmethod
    def _h5_deserialize(cls, h5group, version):
        """
        Read a PPVCube object from a HDF5 file.

        Parameters
        ----------
        h5group: ``h5py.Group``
            Main group to deserialize the PPVCube from.
        version: ``int``
            Version of the PPVCube class to deserialize.

        Returns
        -------
        cube: :class:`~pymses.analysis.ppv.cube.PPVCube`
            new PPVCube instance.
        """
        if version == PPVCube._pymses_h5_version:  # Current PPVCube version
            cam = Camera.from_HDF5(h5group["camera"])
            cube_group = h5group['ppvcube']

            cub = cube_group['data'][...]

            if "name" in cube_group.attrs:
                name = cube_group.attrs['name']
            else:
                name = None

            vel_range = cube_group['velocity_range'][...]
            nvel_bins = cube_group.attrs['nv']

            vel_unit = Unit.from_HDF5(cube_group["vel_unit"])
            Sigma_unit = Unit.from_HDF5(cube_group["sigma_unit"])

        else:
            raise ValueError("Unknown PPVCube version (%d)" % int(version))

        c = cls(cam, sigma_unit=Sigma_unit, vrange=vel_range, nv=nvel_bins, vunit=vel_unit, name=name)
        c._datacube[:] = cub[:]
        return c

__all__ = ["PPVCube"]
