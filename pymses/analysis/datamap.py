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
:mod:`pymses.analysis.datamap` --- Data map module
--------------------------------------------------

"""
import os
import numpy as N

from pymses.utils import HDF5Serializable
from pymses.utils.constants import Unit
from .camera import Camera
from .plot import Plot2D, PlainPNGImage


class _ScalarMapInfo(object):
    def __init__(self, data, unit, camera, data_mask=None, log_scale=None):
        super(_ScalarMapInfo, self).__init__()

        self._camera = camera

        # Mandatory map value unit
        if not isinstance(unit, Unit):
            raise AttributeError("'map_value_unit' is not a valid Unit instance.")
        self._unit = unit

        # Check log_sensitive value
        if log_scale is not None and not isinstance(log_scale, bool):
            raise AttributeError("'log_sensitive' attribute must be a boolean value, if not None.")
        self._is_log = log_scale

        # Check scalar map number of dimensions
        if data.ndim != 2:
            raise AttributeError("Scalar map is not a valid 2D numpy.ndarray instance.")

        self._data = data

        # Save array mask
        if data_mask is not None:
            # Check map mask, if not none, is a valid 2D numpy.ndarray object
            if not isinstance(data_mask, N.ndarray):
                raise AttributeError("'map_mask' is not a valid 2D numpy.ndarray, if not None.")

            # Check mask number of dimensions
            if data_mask.ndim != 2:
                raise AttributeError("Map mask is not a 2D array.")

            # Check min/max mask values
            if data_mask.min() < 0.0 or data_mask.max() > 1.0:
                raise AttributeError("mask array values must lie within [0.0; 1.0].")

            # Mask<->map shape compatibility check
            if self._data.shape != data_mask.shape:
                raise AttributeError("Scalar map shape is not identical to the mask array shape.")

        self._mask = data_mask

    @property
    def is_user_defined_log_scale(self):
        return self._is_log is not None

    @property
    def is_log(self):
        if self._is_log is None:
            return self._camera.log_sensitive

        return self._is_log

    @property
    def data(self):
        return self._data

    @property
    def is_user_defined_mask(self):
        return self._mask is not None

    @property
    def mask(self):
        if self._mask is None:
            mask = self._camera.get_map_mask()

            # Mask<->map shape compatibility check
            if self._data.shape != mask.shape:
                raise ValueError("Scalar map shape is not identical to the mask array shape.")

            return mask

        return self._mask

    @property
    def unit(self):
        return self._unit


class DataMap(HDF5Serializable):
    """
    Data map class

    Parameters
    ----------
    vmap: ``numpy.ndarray``
        raw data of the  default scalar map.
    camera: :class:`~pymses.analysis.camera.Camera`
        camera used to obtain the data
    unit: :class:`~pymses.utils.constants.Unit`
        data value unit.
    mask: ``numpy.ndarray`` or None
        2D mask array to apply to the default scalar map. Default None.
    map_name: ``string``
        name of the map. Default "map".
    """
    _pymses_h5_version = 2

    def __init__(self, vmap, camera, unit, mask=None, map_name="map"):
        super(DataMap, self).__init__()

        # Init map dicts
        self._scalar_maps = {}
        self._points_data = {}
        self._uv_vector_maps = {}

        if not isinstance(camera, Camera):
            raise AttributeError("'camera' attribute is not a valid Camera instance")
        self._camera = camera.copy()

        # Add default map
        self.add_scalar_map(vmap, map_name, unit, map_mask=mask)
        self._default_map_name = map_name

    @property
    def map_size(self):
        return self._camera.get_map_size()

    def add_scalar_map(self, m, map_name, map_value_unit, map_mask=None, log_sensitive=None):
        """
        Add a new scalar 2D map to the DataMap, with a given name.

        Parameters
        ----------
        m: ``numpy.ndarray``
            2D scalar map to add.
        map_name: ``string``
            name of the scalar map.
        map_value_unit: :class: `~pymses.utils.constants.unit.Unit`
            map value unit.
        map_mask: ``numpy.ndarray`` or None
            2D mask array to apply to the added scalar map. Default None, the default value is then taken
            from the Camera.get_map_mask() method.
        log_sensitive: ``bool`` or None.
            Does the datamap need to be displayed in logarithmic scale ? Default None, the default value is then taken
            from the Camera.log_sensitive attribute.
        """
        if not isinstance(map_name, str):
            raise AttributeError("'map_name' attribute must be a valid string.")

        # Checks a scalar map with the same name does not already exist
        if map_name in self._scalar_maps:
            raise AttributeError("A scalar map named '{mname!s}' already exists.".format(mname=map_name))

        # Cast map to a numpy.ndarray
        self._scalar_maps[map_name] = _ScalarMapInfo(N.asarray(m), map_value_unit, self._camera, data_mask=map_mask,
                                                     log_scale=log_sensitive)

    @property
    def nscalar_maps(self):
        return len(self._scalar_maps)

    def _get_scalar_map(self, map_name=None):
        if map_name is not None:
            if not isinstance(map_name, str):
                raise AttributeError("'map_name' attribute is not a valid string.")

            mname = map_name
        else:
            mname = self._default_map_name

        m_info = self._scalar_maps.get(mname, None)
        if m_info is None:
            raise AttributeError("Scalar 2D map with name '{map_name!s}' was not found.".format(map_name=mname))

        return m_info

    @property
    def map(self):
        """
        Backward compatibility : get the default 2D scalar data map array.

        Returns
        -------
        map: 2D ``numpy.ndarray``
            default scalar data map.
        """
        return self.scalar_map()

    def scalar_map(self, map_name=None):
        """
        Fetch scalar map data with a given map name

        Parameters
        ----------
        map_name: ``string`` or None.
            name of the scalar map to fetch. Default None: revert to default scalar map.

        Returns
        -------
        d: 2D ``numpy.ndarray``
            scalar map.
        """
        m_info = self._get_scalar_map(map_name)
        return m_info.data

    def __eq__(self, other):
        if other.nscalar_maps != self.nscalar_maps:
            return False
        for map_name, vmap in self._scalar_maps.items():
            try:
                other_data = other.scalar_map(map_name)
                if not N.allclose(vmap.data, other_data, rtol=1.0e-6):
                    return False
            except AttributeError:
                return False

        # TODO compare point data and vector map data

        if not self._camera != other.camera:
            return False

        return True

    @property
    def camera(self):
        """
        Get the data map camera object

        Returns
        -------
        cam: :class:~pymses.analysis.camera.Camera`
            data map camera instance
        """
        return self._camera

    @property
    def map_unit(self):
        """
        Backward compatibility : get the default 2D scalar data map value unit.

        Returns
        -------
        u: :class:~pymses.utils.constants.unit.Unit`
            default scalar data map value unit.
        """
        return self._get_scalar_map().unit

    def value_range(self, map_name=None, vrange=None, fraction=1.0, log_scale=None, verbose=False):
        """
        Map range computation function. Computes the linear/log scale map value range.

            * if a user-defined vrange is given, then it is used to compute the map range values, e.g. (vmin, None), or
              (vmin, vmax), or (None, vmax).
            * if not, the map range values is computed from a fraction (percent) of the total value
              of the map parameter. the min. map range value is defined as the value below which there
              is a fraction of the map (default 1 %)

        Parameters
        ----------
        map_name   : name of the scalar map from which the map range values are computed.
        vrange     : user-defined map value range.
        fraction   : fraction of the total map values below the min. map range (in percent). Default 1 %.
        log_scale  : map values displayed in log-scale ? Default None.

        Returns
        -------
        map_range : (``float``, ``float``) ``tuple``
            the map range values (vmin, vmax)
        """
        # Get a copy of the transposed 2D scalar map
        m_info = self._get_scalar_map(map_name)
        # mt = m_info.data.copy().transpose()
        mt = m_info.data
        mask = m_info.mask

        is_log = m_info.is_log
        if log_scale is not None and isinstance(log_scale, bool):
            is_log = log_scale

        if mask.max() == 0.0:  # Map mask is completely null => region outside of simulation domain
            return None, None

        vmin = None
        vmax = None

        if vrange is not None:  # Unfold vrange into a (vmin, vmax) tuple
            if not isinstance(vrange, tuple):
                raise AttributeError("vrange parameter must be a (vmin, vmax) tuple. Got %s" % type(vrange))
            if vrange[1] is not None:
                vmax = vrange[1]

            if vrange[0] is not None:
                vmin = vrange[0]

        if vmin is None or vmax is None:
            if mask.min() > 0.0:  # Map mask has no transparent pixel => get all values
                range_mvalues = mt.reshape(-1)[:]
            else:
                range_mvalues = mt[mask > 0.0]

            if vmax is None:  # Default vmax = map maximum value
                vmax = range_mvalues.max()

            # Fallback values to (None, None) when all map values are negative or null in log scale view
            if is_log and vmax <= 0.0:
                raise ValueError("Cannot compute map range when all log-scaled map values are negative.")

            if vmin is None:
                if fraction < 0.0 or fraction > 100.0:
                    raise AttributeError("fraction parameter must be in the range [0.0, 100.0]. "
                                         "Got {frac:g}".format(frac=fraction))
                frac = fraction / 100.0

                if frac == 0.0:
                    if is_log:
                        vmin = range_mvalues[range_mvalues > 0.0].min()
                    else:
                        vmin = range_mvalues.min()
                elif frac == 1.0:
                    vmin = vmax
                else:
                    if is_log:  # Sort strictly positive values of the map value list
                        values = N.sort(range_mvalues[range_mvalues > 0.0])
                    else:  # Sort values of the masked map
                        values = N.sort(range_mvalues)
                    weights = values.copy()

                    wrange = weights[-1] - weights[0]

                    # Flat map case
                    if wrange == 0.0:
                        vmin = values[0]
                    else:
                        weights = (weights - weights[0]) / wrange
                        # vmin = N.percentile((weights-weights[0])/wrange, q=frac*100.0, interpolation='higher')
                        cumval = N.cumsum(weights)
                        cumval /= cumval[-1]
                        mask = (cumval >= frac)
                        vmin = values[mask][0]

        if vmin >= vmax:
            raise AttributeError("min. value ({min:g}) is greater or equal to max. value ({max:g}).".format(min=vmin,
                                                                                                            max=vmax))

        if verbose:
            if is_log:
                scaling = "logarithmic"
            else:
                scaling = "linear"
            if map_name is None:
                mname = self._default_map_name
            else:
                mname = map_name
            print("Map ('{map_name!s}') value range is : [{min:g}, {max:g}] ({sc!s} scale)".format(map_name=mname,
                                                                                                   min=vmin, max=vmax,
                                                                                                   sc=scaling))
        return (vmin, vmax), is_log

    def save_PNG(self, img_gen=None, map_name=None, img_fname=None, vrange=None, fraction=1.0, verbose=False):
        """
        Convert the map into a PIL Image and save it into a PNG file.

        Parameters
        ----------
        img_gen: TODO
            TODO
        map_name: ``string`` or None
            name of the scalar map (among scalar map dictionary values) to save as a PNG file. Default None: revert to
            default scalar map.
        img_fname: ``string``
            PNG image file name or path
        vrange: ``tuple`` or `None`
            linear-scale map value range (vmin, vmax) tuple used for clipping before saving the PNG file. Default is
            None. When left to None, the map value range is computed automatically. One boundary may be set by
            defining (vmin, None) or (None, vmax).
        fraction: ``float``
            Used by automatic map value range computation. Fraction of the cumulated map values (linear or log scaled)
            below the min. map range (in percent). Default : 1 %.
        verbose: ``bool``
            Verbose processing ? Default : False.
        """
        # Get or init. a PNG image generator
        image_generator = img_gen
        if img_gen is None:
            image_generator = PlainPNGImage()

        # Get the transposed 2D scalar map
        m_info = self._get_scalar_map(map_name)
        m = m_info.data.transpose()

        vrange_out, is_log = self.value_range(map_name, vrange, fraction, image_generator.force_log_scale, verbose)

        return image_generator(m, m_info.mask.transpose(), vrange_out, is_log, img_fname, verbose)

    def save_plot(self, plot_gen=None, map_name=None, plot_fname=None, vrange=None, fraction=1.0, verbose=False):
        """
        Convert the map into a matplotlib plot and save it into a PNG file, if required.

        Parameters
        ----------
        plot_gen: matplotlib Figure generator
        map_name: TODO
        plot_fname: ``string`` or None
            PNG image file name or path. Default None : do not save the plot into a PNG file.
        vrange: ``tuple`` or `None`
            linear-scale map value range (vmin, vmax) tuple used for clipping before saving the PNG file. Default is
            None. When left to None, the map value range is computed automatically. One boundary may be set by
            defining (vmin, None) or (None, vmax).The range values are expressed in 'map_unit', if defined. If
            'map_unit' is not defined, the range value is expressed in base map value unit.
        fraction: ``float``
            Used by automatic map value range computation. Fraction of the cumulated map values (linear or log scaled)
            below the min. map range (in percent). Default : 1 %.
        verbose: ``bool``
            Verbose processing ? Default : False.
        """
        # Get or init. a plot generator
        plot_generator = plot_gen
        if plot_gen is None:
            plot_generator = Plot2D()

        # Get a the transposed 2D scalar map
        m_info = self._get_scalar_map(map_name)
        m = m_info.data.transpose()

        vrange_out, is_log = self.value_range(map_name, vrange, fraction, plot_generator.force_log_scale, verbose)

        # Get u/v axes info
        axes_info = self._camera.get_uvaxes_edges_labels(axis_unit=plot_generator.axis_unit,
                                                         force_zero_centered_axis_values=plot_generator.force_zero_axis)

        return plot_generator(m, m_info.unit, axes_info, vrange_out, is_log, plot_fname, verbose)

    def save_FITS(self, fits_fname, map_name=None,  axis_unit=None, map_unit=None):
        r"""
        Function that saves the DataMap into a FITS file.

        Parameters
        ----------
        fits_fname : the FITS file path in which the DataMap is to be saved
        axis_unit: :class:`~pymses.utils.constants.unit.Unit`
            U/V-axis output unit
        map_unit: :class:`~pymses.utils.constants.unit.Unit`
            map value output unit
        """
        try:
            from astropy.io import fits
            from astropy import wcs
        except ImportError:
            raise ImportError("astropy package is not available. It is mandatory to save FITS files.")

        hdr = self._camera.get_fits_header(axis_unit=axis_unit)

        # Get the transposed 2D scalar map
        m_info = self._get_scalar_map(map_name)
        vmap = m_info.data.transpose().copy()
        if map_unit is not None:
            hdr["BUNIT"] = map_unit.name
            vmap *= self._unit.express(map_unit)
        else:
            hdr["BUNIT"] = self._unit.name

        fits_dir = os.path.dirname(fits_fname)
        print("Saving image into FITS file '%s'" % fits_fname)
        if not os.path.isdir(fits_dir):
            os.makedirs(fits_dir)
        hdu = fits.PrimaryHDU(vmap, header=hdr)
        hdulist = fits.HDUList([hdu])

        # Write DataMap to FITS file, overwriting it if the file already exists
        hdulist.writeto(fits_fname, clobber=True)
        print("FITS File '%s' saved." % fits_fname)

    def _h5_serialize(self, h5group, **kwargs):
        """
        Serialize the DataMap object into a HDF5 file.

        Parameters
        ----------
        h5group: ``h5py.Group``
            Main group to serialize the DataMap into
        float32: ``bool``
            save map data with dtype=float32 instead of float64 ? default False
        save_mask: ``bool``
        """
        float32 = kwargs.get("float32", False)

        # Save camera object
        cam_group = h5group.create_group("camera")
        self._camera.save_HDF5(cam_group)

        # Save all scalar map data
        scalar_group = h5group.create_group("scalar_maps")
        scalar_group.attrs['nmaps'] = len(self._scalar_maps)
        imap = 0
        for map_name, smi in self._scalar_maps.items():
            # Create map group
            map_group = scalar_group.create_group("map_{map_number:d}".format(map_number=imap))

            # Save map name
            map_group.attrs['map_name'] = str(map_name)
            if map_name == self._default_map_name:  # Store default map index in parent group attributes
                scalar_group.attrs['default_map_index'] = imap

            # Save map data array
            if float32:
                # use float32 type to save disk space :
                smap = smi.data.astype("float32")
            else:
                smap = smi.data
            map_group.create_dataset("data", data=smap, compression='gzip', compression_opts=9)  # Best compression
            # Save optional (user-defined) log-scale attribute
            if smi.is_user_defined_log_scale:
                map_group.attrs['is_log'] = smi.is_log

            # Save map mask array, only if user-defined
            if smi.is_user_defined_mask:
                if float32:
                    smask = smi.mask.astype("float32")
                else:
                    smask = smi.mask
                map_group.create_dataset("mask", data=smask, compression='gzip', compression_opts=9)  # Best compression

            # Save map meta-info : value range
            map_range = N.array([N.min(smap), N.max(smap)])
            map_group.create_dataset("value_range", data=map_range)

            # Save scalar map value unit
            unit_group = map_group.create_group("unit")
            smi.unit.save_HDF5(unit_group)

            imap += 1

        # Save all point data
        if len(self._points_data) > 0:
            point_data_group = h5group.create_group("point_data")
            point_data_group.attrs['ndatasets'] = len(self._points_data)
            # TODO : save points data

        # Save all vector map data
        if len(self._uv_vector_maps) > 0:
            vector_group = h5group.create_group("vector_maps")
            vector_group.attrs['nmaps'] = len(self._uv_vector_maps)
            # TODO : save (u,v) vector map data

    @classmethod
    def _h5_deserialize(cls, h5group, version):
        """
        Read a DataMap object from a HDF5 file.

        Parameters
        ----------
        h5group: ``h5py.Group``
            Main group to deserialize the DataMap from.
        version: ``int``
            Version of the DataMap class to deserialize.

        Returns
        -------
        dm: :class:`~pymses.analysis.datamap.DataMap`
            new Datamap instance
        """
        if version == 1:
            cam = Camera.from_HDF5(h5group["camera"])
            map_group = h5group['map']

            map = map_group['data'][...]

            if 'name' in map_group.attrs:
                name = map_group.attrs['name']
            else:
                name = "map"

            # mask = map_group['mask'][...]
            # map_size = map_group['size'][...]
            # map_range = map_group['value_range'][...]

            value_unit = Unit.from_HDF5(map_group['unit'])

            m = cls(map, cam, value_unit, map_name=name)
            return m

        elif version == DataMap._pymses_h5_version:  # Current DataMap version (2)
            # Read camera
            cam = Camera.from_HDF5(h5group["camera"])

            # Read all scalar map data
            scalar_group = h5group["scalar_maps"]
            nmaps = scalar_group.attrs['nmaps']
            default_map_index = scalar_group.attrs['default_map_index']

            # Read default scalar map
            map_group = scalar_group["map_{map_number:d}".format(map_number=default_map_index)]
            mname = map_group.attrs['map_name']
            scal_map = map_group['data'][...]
            if 'mask' in map_group:
                map_mask = map_group['mask'][...]
            else:
                map_mask = None
            value_unit = Unit.from_HDF5(map_group['unit'])
            dm = cls(scal_map, cam, value_unit, mask=map_mask, map_name=mname)

            for imap in range(nmaps):  # Read all additional scalar maps
                if imap == default_map_index:
                    continue
                map_group = scalar_group["map_{map_number:d}".format(map_number=imap)]
                mname = map_group['map_name']

                scal_map = map_group['data'][...]
                if 'mask' in map_group:
                    map_mask = map_group['mask'][...]
                else:
                    map_mask = None
                value_unit = Unit.from_HDF5(map_group['unit'])
                if 'is_log' in map_group.attrs:
                    log_scale = map_group.attrs['is_log']
                else:
                    log_scale = None

                dm.add_scalar_map(scal_map, mname, value_unit, map_mask, log_sensitive=log_scale)

            # TODO : read point data

            # TODO : read (u,v) vector map data
            return dm
        else:
            raise ValueError("Unknown DataMap version (%d)" % int(version))

    @classmethod
    def load_legacy_HDF5_datamap(cls, h5fname):
        """
        Read a DataMap object from a legacy (PyMSES v<=4.0.0) HDF5 file.

        Parameters
        ----------
        h5fname: ``string``
            HDF5 file name.

        Returns
        -------
        dmap: :class:`~pymses.analysis.datamap.DataMap`
            DataMap instance
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
        try:
            h5f = T.File(h5fname, 'r')
            h5g = h5f.get_node("/")

            # Old datamap version (PyTables syntax)
            map_name = h5g.name.read()
            cam = Camera._load_legacy_cam(h5f.get_node("/camera"))
            map_group = h5g.map
            map = map_group.map.read()
            # Mandatory map value unit
            value_unit = Unit._load_legacy_HDF5_unit(map_group.unit)
            # Mandatory size unit unit => save it into Camera instance
            cam.size_unit = Unit._load_legacy_HDF5_unit(map_group.length_unit)

            # vrange_path = map_group._g_join("value_range")
            # if h5f.__contains__(vrange_path):
            #     map_range = map_group.value_range.read()
            # else:
            #     map_range = None
            #
            # msize_path = map_group._g_join("size")
            # if h5f.__contains__(msize_path):
            #     map_size = map_group.size.read()
            # else:
            #     map_size = None

            dmap = cls(map, cam, value_unit, map_name=map_name)
        except Exception as exc:
            raise IOError("HDF5 I/O error : %s" % exc)
        finally:
            h5f.close()

        return dmap

    @classmethod
    def value_range_from_multiple_HDF5(cls, h5fname_iter, fraction=None):
        """
        Get the common value range from a DataMap HDF5 file iterator

        Parameters
        ----------
        h5fname_iter: DataMap HDF5 filename iterator

        Returns
        -------
        vrange: ``tuple`` of ``float``
            (vmin, vmax) global value range
        """
        all_vmin = None
        all_vmax = None
        for h5fname in h5fname_iter:
            d = DataMap.from_HDF5(h5fname)

            if fraction is not None:
                vrange, is_log = d.value_range(fraction=fraction)
            else:
                vrange, is_log = d.value_range()
            vmin, vmax = vrange

            if all_vmin is None:
                all_vmin = vmin
            elif vmin < all_vmin:
                all_vmin = vmin
            if all_vmax is None:
                all_vmax = vmax
            elif vmax > all_vmax:
                all_vmax = vmax

        return all_vmin, all_vmax


__all__ = ["DataMap"]
