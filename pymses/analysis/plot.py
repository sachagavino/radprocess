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
:mod:`pymses.analysis.plot` --- plot module
-------------------------------------------

"""
import os
from PIL import Image as I
import numpy as N
from matplotlib import pyplot as P
from matplotlib.figure import Figure
from matplotlib.ticker import IndexLocator, FormatStrFormatter


from pymses.utils.colors import Colormaps
from pymses.utils.filename import FileUtil
from pymses.utils.constants import Unit


def _apply_log_scale(vmap, vrange, verbose=False):
    """
    Apply log-scaling to a data array

    Parameters
    ----------
    vmap: (nx, ny) ``numpy.ndarray``
        map values 2D array
    verbose: ``bool``
        is verbose ? default False.
    """
    # Apply log-scale map
    neg_null_mask = vmap <= 0.0

    if neg_null_mask.all():  # No positive values ?
        raise ValueError("Cannot apply log scale : all map values are <= 0 !")

    if neg_null_mask.any():  # Clip negative or null map values to apply log scale safely
        if vrange[0] is not None:
            non_null_min_map = 10.0 ** vrange[0]  # value range minimum log value
        else:  # Should never happen -> all map values must be null if mask is null
            non_null_min_map = N.min(vmap[~neg_null_mask]) / 10.0  # 1/10th of the lowest positive map value

        # Clipping
        vmap[neg_null_mask] = non_null_min_map
        if verbose:
            print("Warning: %d values <= 0 were replaced by min_map/10 = %g in order to apply log scale" % \
                    (N.sum(neg_null_mask), non_null_min_map))

    # Inplace log-scaling
    N.log10(vmap, out=vmap)


class AbstractDisplayGenerator(object):
    """
    2D data view generator abstract class
    """
    def __init__(self, cmap="Viridis", discrete=False, force_log_scale=None):
        self._force_log_scale = None
        self.force_log_scale = force_log_scale
        self._discrete = discrete
        self._cmap = cmap

    @property
    def is_discrete_colormap(self):
        return self._discrete

    @property
    def force_log_scale(self):
        return self._force_log_scale

    @force_log_scale.setter
    def force_log_scale(self, force_log_scale):
        if force_log_scale is not None and not isinstance(force_log_scale, bool):
            raise AttributeError("'force_log_scale' is not a valid boolean value, if not None.")
        self._force_log_scale = force_log_scale

    def _get_colormap(self, vrange_out, is_log_scale, verbose=False):
        """
        Compute colormap
        """
        vmin, vmax = vrange_out

        if is_log_scale:
            vmin = N.log10(vmin)
            vmax = N.log10(vmax)

        if self._discrete:
            self._colormap, vrange = Colormaps.get_cmap(self._cmap, discrete=True, value_range=(vmin, vmax))
            vmin, vmax = vrange
        else:
            self._colormap = Colormaps.get_cmap(self._cmap)

        vrange_cmap = (vmin, vmax)

        if verbose:
            print("Colormap range is : [%g %g]" % (vmin, vmax))

        return vrange_cmap


class PlainPNGImage(AbstractDisplayGenerator):
    """
    Pillow 2D image generator

    Parameters
    ----------
    cmap: ``string`` or ``matplotlib.colors.Colormap``
        applied colormap to generate the PIL Image.
    discrete: ``bool``
        Show only discrete values in colormap. Default : False.
    force_log_scale: ``bool`` or `None`
        True if the image file is a log-scale view of a 2D data map, otherwise False. default None. If None, the
        `log_sensitive` parameter of the Datamap :class:`~pymses.analysis.camera.Camera` instance is used.
    alpha_mask: ``bool``
        Use the map mask to generate the PNG file alpha band ? Default : True.
    """
    def __init__(self, cmap="Viridis", discrete=False, force_log_scale=None, alpha_mask=True):
        super(PlainPNGImage, self).__init__(cmap, discrete, force_log_scale)
        self._alpha_mask = alpha_mask

        self._nx = -1
        self._ny = -1
        self._nx_mask = -1
        self._ny_mask = -1
        self._R_band = None
        self._G_band = None
        self._B_band = None
        self._A_band = None

    def _put_pil_img(self, map_01, mask_01=None):
        """
        Put a 2D [0,1] map (+ optional) mask array into a PIL image RGB(+A) bands
        :param map_01:
        :param mask_01:
        :return:
        """
        ny, nx = map_01.shape

        change_size = False
        if self._nx != nx or self._ny != ny:
            change_size = True
            self._nx = nx
            self._ny = ny

        m = N.clip((self._colormap(map_01)*255.0).astype('i'), 0, 255).reshape(nx * ny, 4)

        if self._R_band is None or change_size:
            self._R_band = I.new("L", (nx, ny))
        self._R_band.putdata(m[:, 0])
        if self._G_band is None or change_size:
            self._G_band = I.new("L", (nx, ny))
        self._G_band.putdata(m[:, 1])
        if self._B_band is None or change_size:
            self._B_band = I.new("L", (nx, ny))
        self._B_band.putdata(m[:, 2])

        if mask_01 is not None and self._alpha_mask:
            change_msk_size = False
            if self._nx_mask != nx or self._ny_mask != ny:
                change_msk_size = True
                self._nx_mask = nx
                self._ny_mask = ny
            map_mask = N.clip((256.0 * mask_01).astype('i'), 0, 255)
            map_mask = map_mask.reshape(nx * ny)
            if self._A_band is None or change_msk_size:
                self._A_band = I.new("L", (nx, ny))
            self._A_band.putdata(map_mask)
            out_img = I.merge("RGBA", (self._R_band, self._G_band, self._B_band, self._A_band))
        else:
            out_img = I.merge("RGB", (self._R_band, self._G_band, self._B_band))

        # Flip top <-> bottom (image origin corner is top left corner)
        out_img = out_img.transpose(I.FLIP_TOP_BOTTOM)

        return out_img

    def __call__(self, value_map, mask, value_range, is_log_scale, img_fname=None, verbose=False):
        """
        Convert the map into a PIL Image and save it into a PNG file.

        Parameters
        ----------
        value_map: TODO
        mask: TODO
        value_range: ``tuple`` or `None`
            linear-scale map value range (vmin, vmax) tuple used for clipping before saving the PNG file. Default is
            None. When left to None, the map value range is computed automatically. One boundary may be set by
            defining (vmin, None) or (None, vmax).
        is_log_scale: ``bool``
            True if the image file is a log-scale view of the map, otherwise False.
        img_fname: ``string``
            PNG image path. Default None: do no save the Pillow Image as a PNG file.
        verbose: ``bool``
            Verbose processing ? Default : False.

        Returns
        -------
        out_img: Pillow image.
        """
        vrange_cmap = self._get_colormap(value_range, is_log_scale, verbose)
        vmap = value_map.copy()

        # Log scale ?
        if is_log_scale:
            _apply_log_scale(vmap, vrange_cmap, verbose)

        # Clipping
        N.clip(vmap, vrange_cmap[0], vrange_cmap[1], out=vmap)

        vmap = (vmap - vrange_cmap[0]) / (vrange_cmap[1] - vrange_cmap[0])
        if mask.min() == 1.0:
            map_mask = None
        else:
            map_mask = mask.copy()

        out_img = self._put_pil_img(vmap, mask_01=map_mask)

        if img_fname is not None:
            # Save image
            path = FileUtil.new_filepath(img_fname, append_extension=FileUtil.PNG_FILE)
            print("Saving img into '%s'" % path)
            out_img.save(path)

        return out_img


class Plot2DFigure(Figure):
    def remove(self):
        pass

    def display_vector_field(self, vecx, vecy):
        pass


class Plot2D(AbstractDisplayGenerator):
    """
    Matplotlib 2D plot generator

    Parameters
    ----------
    cmap: ``string`` or ``matplotlib.colors.Colormap``
        applied colormap to generate the PIL Image.
    discrete: ``bool``
        Show only discrete values in colormap. Default : False.
    force_log_scale: ``bool`` or `None`
        True if the image file is a log-scale view of a 2D data map, otherwise False. default None. If None, the
        `log_sensitive` parameter of the Datamap :class:`~pymses.analysis.camera.Camera` instance is used.
    axis_unit: :class:`~pymses.utils.constants.unit.Unit`
        map axis size unit used in output plot.
    map_unit: :class:`~pymses.utils.constants.unit.Unit`
        map value unit used in output plot.
    force_zero_centered_axis_values: ``bool``
        Attribute forwarded to :meth:`Camera.get_uvaxes_edges_labels()` . Default False.
    """
    def __init__(self, cmap="Viridis", discrete=False, force_log_scale=None, axis_unit=None, map_unit=None,
                 force_zero_centered_axis_values=False):
        super(Plot2D, self).__init__(cmap, discrete, force_log_scale)

        self._axis_unit = None
        self.axis_unit = axis_unit

        self._map_value_unit = None
        self.map_value_unit = map_unit
        self._map_extent = None
        self._current_aximage = None # Current matplotlib.image.AxisImage instance

        self._force_zerocen_axvalues = False
        self.force_zero_axis = force_zero_centered_axis_values

    @property
    def axis_unit(self):
        return self._axis_unit

    @axis_unit.setter
    def axis_unit(self, new_axis_unit):
        if new_axis_unit is not None and not isinstance(new_axis_unit, Unit):
            raise AttributeError("'axis_unit' is not a valid Unit instance.")
        self._axis_unit = new_axis_unit

    @property
    def map_value_unit(self):
        return self._map_value_unit

    @map_value_unit.setter
    def map_value_unit(self, new_map_value_unit):
        if new_map_value_unit is not None and not isinstance(new_map_value_unit, Unit):
            raise AttributeError("'map_value_unit' is not a valid Unit instance.")
        self._map_value_unit = new_map_value_unit

    @property
    def force_zero_axis(self):
        return self._force_zerocen_axvalues

    @force_zero_axis.setter
    def force_zero_axis(self, new_force_zero_center):
        if not isinstance(new_force_zero_center, bool):
            raise AttributeError("'force_zero_axis' value must be boolean.")
        self._force_zerocen_axvalues = new_force_zero_center

    @property
    def current_axis_image(self):
        return self._current_aximage

    @property
    def image_extent(self):
        return self._map_extent.copy()

    def __call__(self, value_map, value_unit, axes_info, value_range, is_log_scale, plot_fname=None, verbose=False):
        """
        Plot map using matplotlib.pyplot.imshow()

        Parameters
        ----------
        :param value_map:
        :param value_unit:
        :param axes_info:
        :param value_range:
        :param is_log_scale:
        :param plot_fname:
        :param verbose:

        Returns
        -------
        """
        vrange_cmap = self._get_colormap(value_range, is_log_scale, verbose)
        vmap = value_map.copy()

        # Log scale ?
        if is_log_scale:
            _apply_log_scale(vmap, vrange_cmap, verbose)

        # Clipping
        N.clip(vmap, vrange_cmap[0], vrange_cmap[1], out=vmap)

        # Map value unit conversion
        if self._map_value_unit is not None and isinstance(self._map_value_unit, Unit):
            map_value_factor = value_unit.express(self._map_value_unit)
            vmap *= map_value_factor
            map_unit_label = self._map_value_unit.latex
        else:
            map_unit_label = value_unit.latex

        # Get u/v axes labels and pixel edge coordinates
        uinfo, vinfo = axes_info
        _u_axisname, _u_axisunit, u_label_latex, uedges = uinfo
        _v_axisname, _v_axisunit, v_label_latex, vedges = vinfo

        # Plot matplotlib.pyplot Figure
        P.figure(FigureClass=Plot2DFigure)
        self._map_extent = [uedges[0], uedges[-1], vedges[0], vedges[-1]]
        self._current_aximage = P.imshow(vmap, cmap=self._colormap, vmin=vrange_cmap[0], vmax=vrange_cmap[1],
                                         origin='lower', extent=self._map_extent, interpolation='none')

        # Set axis labels
        P.xlabel(u_label_latex)
        P.ylabel(v_label_latex)

        # Pretty user-defined colorbar
        if is_log_scale:
            fo = FormatStrFormatter("$10^{%d}$")
            offset = N.ceil(vrange_cmap[0]) - vrange_cmap[0]
            lo = IndexLocator(1.0, offset)
            cb = P.colorbar(ticks=lo, format=fo)
        else:
            if self._discrete:
                fo = FormatStrFormatter("%d")
                ncol = int(N.round(vrange_cmap[1] - vrange_cmap[0]))
                ti = N.linspace(vrange_cmap[0] + 0.5, vrange_cmap[1] - 0.5, ncol)
                cb = P.colorbar(format=fo, ticks=ti)
            else:
                cb = P.colorbar()

        # Set colorbar label
        cb.set_label(map_unit_label)

        fig = P.gcf()
        # Automatically adjust layout to fit the figure canvas
        fig.tight_layout()

        # Save image
        if plot_fname is not None:
            fname = FileUtil.new_filepath(plot_fname, append_extension=FileUtil.PNG_FILE)
            print("Saving plot into '%s'" % fname)
            fig.savefig(fname)

        return fig


__all__ = ['Plot2D', 'PlainPNGImage']
