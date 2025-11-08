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
:mod:`pymses.utils.colormaps` --- Custom PyMSES colormaps module
----------------------------------------------------------------

"""
import numpy as N
from matplotlib import pyplot as P
from matplotlib import cm
from matplotlib.colors import Colormap, LinearSegmentedColormap
from .custom import custom_cmaps
from .pus import PUS_cmaps


class Colormaps(object):
    _custom_maps = {"BlackBlueWhiteRed": custom_cmaps['BlackBlueWhiteRed'],
                    "BlackBlueYellowRed": custom_cmaps['BlackBlueYellowRed'],
                    "BlueYellowTanager": custom_cmaps['BlueYellowTanager'],
                    "BlackPurpWhiBlGreen": custom_cmaps['BlackPurpWhiBlGreen'],
                    "Dark_green": custom_cmaps['Dark_green'],
                    "JetBlack": custom_cmaps['JetBlack'],
                    "Viridis": PUS_cmaps['Viridis'],
                    "Inferno": PUS_cmaps['Inferno'],
                    "Plasma": PUS_cmaps['Plasma'],
                    "Magma": PUS_cmaps['Magma'],
                    "Viridis_dark": custom_cmaps['Viridis_dark'],
                    "Dark_thermal": custom_cmaps['Dark_thermal'],
                    "WhiteBlueYellow": custom_cmaps['WhiteBlueYellow']
    }

    @classmethod
    def cmap_list(cls):
        l = list(Colormaps._custom_maps.keys()) + list(cm.cmap_d.keys())
        l.sort()
        return l

    @classmethod
    def get_cmap(cls, cmap, discrete=False, value_range=None):
        """
        Get a matplotlib.cm.colors.Colormap instance from the PyMSES custom colormap catalog or the default matplotlib
        catalog.


        Parameters
        ----------
        cmap: ``string`` or :class:`matplotlib.cm.colors.Colormap` instance.
            name of the colormap or Colormap instance.
        discrete: ``bool``
            Is the colormap needs to be made of discrete values ? Use with 'value_range' attribute value. Default false.
        value_range: ``tuple`` of ``float``
            User-defined (vmin, vmax) value range tuple. For use with discrete=True. Default None.

        Examples
        --------
        TODO

        >>> cmap = colormaps.get_cmap('Viridis')

        Returns
        -------
        colormap: :class:`matplotlib.cm.colors.Colormap` instance
            required colormap
        value_range: ``tuple`` of ``float``
            If discrete colormap was requested, a new (vmin_new, vmax_new) value range tuple is also returned.
        """
        if isinstance(cmap, Colormap):
            colormap = cmap
        else:
            # Sanity check
            if not isinstance(cmap, str):
                raise AttributeError("'name' must be a valid colormap name (string) !")

            if cmap in Colormaps._custom_maps:
                # Try to find colormap in PyMSES custom catalog
                colormap = Colormaps._custom_maps[cmap]
            else:
                # Try to find the colormap in matplotlib colormap catalog
                try:
                    colormap = cm.get_cmap(cmap)
                except ValueError:
                    raise AttributeError(
                        "'%s' colormap neither found in custom PyMSES colormap table :\n%s\n nor in "
                        "matplotlib colormaps :\n%s" % (cmap, list(Colormaps._custom_maps.keys()), list(cm.cmap_d.keys())))
        if discrete and value_range is not None:
            # Modifies the linear segmented colormap in order to show only a restricted number of colors for the
            # required value range
            vmin, vmax = value_range

            ncols = int(round(vmax - vmin)) + 1
            colormap = LinearSegmentedColormap('%s_discrete' % colormap.name, colormap._segmentdata, ncols)
            vmin_out = vmin - 0.5
            vmax_out = vmax + 0.5
            new_vrange = (vmin_out, vmax_out)

            return colormap, new_vrange
        else:
            return colormap

    @classmethod
    def show_colormaps(cls, cm_list=None, sort=True):
        cml = cls.cmap_list()
        if cm_list is not None:
            cml = [cm for cm in cm_list if cm in cml]
        if sort:
            cml.sort()

        nrows = len(cml)
        gradient = N.linspace(0, 1, 256)
        gradient = N.vstack((gradient, gradient))

        fig, axes = P.subplots(nrows=nrows)
        fig.subplots_adjust(top=0.95, bottom=0.01, left=0.2, right=0.99)
        axes[0].set_title("Colormaps", fontsize=14)

        for ax, name in zip(axes, cml):
            ax.imshow(gradient, aspect='auto', cmap=cls.get_cmap(name))
            pos = list(ax.get_position().bounds)
            x_text = pos[0] - 0.01
            y_text = pos[1] + pos[3] / 2.
            fig.text(x_text, y_text, name, va='center', ha='right', fontsize=10)

        # Turn off *all* ticks & spines, not just the ones with colormaps.
        for ax in axes:
            ax.set_axis_off()

        P.show()


__all__ = ["Colormaps"]
