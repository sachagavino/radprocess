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
import numpy as N
from matplotlib.cm import get_cmap


class ColorLinesTransferFunction(object):
    def __init__(self, vbounds, cmap="hsv", alpha=None, gamma=0.6):
        r"""A complete set of transfer functions for standard color-mapping.

        This is the best and easiest way to set up volume rendering.  It
        creates field tables for all three colors, their alphas, and has
        support for sampling color maps and adding independent color values at

        Parameters
        ----------
        vbounds : tuple of floats
            The min and max for the transfer function.	Values below or above
            these values are discarded.
        cmap : string
            Matplotlib Colormap key.
        alpha : func
            alpha function (default None)
        gamma: float
            gamma value (emissivity)

        """
        self.vmin, self.vmax = vbounds
        assert (self.vmin < self.vmax)
        self.cmap = get_cmap(cmap)

        # Opacity function
        if alpha is None:
            def alph_func(v):
                if v <= self.vmin:
                    return 1.0E-3
                elif v >= self.vmax:
                    return 1.0
                else:
                    return 10. ** (-3.0 + (v - self.vmin) * 3.0 / (self.vmax - self.vmin))

            self.alpha = alph_func
        else:
            self.alpha = alpha

        self.gamma = gamma

        self.vlines = []
        self.wlines = []
        self.red_lines = []
        self.green_lines = []
        self.blue_lines = []
        self.alpha_lines = []

    def similar(self, tf):
        """
        Draftly test if a transfer functions is roughly equal to an other one, just to know in the rt processing code
        if we need to reload the transfer functions or not.
        """
        return self.gamma == tf.gamma and self.vlines == tf.vlines and self.wlines == tf.wlines

    def add_line(self, value, width=0.1):
        r"""Add a Gaussian distribution to the color lines transfer function.

        Parameters
        ----------
        value : float
            The center of the gaussian distribution
        width : float
            The gaussian width

        height : list of 4 float
            The peak height (:math:`h` in the above equation.)	Note that while
            values greater 1.0 will be accepted, the values of the transmission
            function are clipped at 1.0.  This must be a list, and it is in the
            order of (red, green, blue, alpha).
        height : float
            The peak height (:math:`h` in the above equation.)	Note that while
            values greater 1.0 will be accepted, the values of the transmission
            function are clipped at 1.0.

        Examples
        --------
        This adds a red spike.

        >>> cltf = ColorLinesTransferFunction( (-5.0, 2.0) )
        >>> cltf.add_line(-2.0, 0.2)
        """
        if ((value < self.vmin) or (value > self.vmax)):
            print("Wrong line position, add_line() ignored")
        else:
            for v, w in zip(self.vlines, self.wlines):
                assert abs(value - v) > (3. * (w + width)), \
                    "Error, intersecting gausian distribution : (x0, sigma0) = (%f, %f) and (x1, sigma1) = (%f, %f)" % (
                        value, width, v, w)
            self.vlines.append(value * 1.)
            self.wlines.append(width * 1.)
            self.alpha_lines.append(self.alpha(value))

            x = (value - self.vmin) / (self.vmax - self.vmin)
            r, g, b, a = self.cmap(x)
            self.red_lines.append(r)
            self.green_lines.append(g)
            self.blue_lines.append(b)

    @property
    def vwrgba(self):
        v = N.asarray(self.vlines)
        ind = N.argsort(v)
        vs = v[ind]
        ws = N.asarray(self.wlines)[ind]
        rs = N.asarray(self.red_lines)[ind]
        gs = N.asarray(self.green_lines)[ind]
        bs = N.asarray(self.blue_lines)[ind]
        als = N.asarray(self.alpha_lines)[ind]
        return (vs, ws, rs, gs, bs, als)


# def plot(self, fname):
#		 r"""Save an image file of the color transfer function.
#
#		 Parameters
#		 ----------
#		 fname : string
#			 The filename where to save the ColorLinesTransferFunction plot
#
#		 Examples
#		 --------
#
#		 >>> tf = ColorLinesTransferFunction( (-10.0, -5.0) )
#		 >>> tf.add_line()
#		 >>> tf.plot("cmap.png")
#		 """
#		 from matplotlib import pyplot
#		 from matplotlib.ticker import FuncFormatter
#		 
#		self.x = N.linspace(x_bounds[0], x_bounds[1], nbins).astype('float64')
#		 
#		pyplot.clf()
#		 ax = pyplot.axes()
#		 i_data = N.zeros((self.alpha.x.size, self.funcs[0].y.size, 3))
#		 i_data[:,:,0] = N.outer(N.ones(self.alpha.x.size), self.funcs[0].y)
#		 i_data[:,:,1] = N.outer(N.ones(self.alpha.x.size), self.funcs[1].y)
#		 i_data[:,:,2] = N.outer(N.ones(self.alpha.x.size), self.funcs[2].y)
#		 ax.imshow(i_data, origin='lower')
#		 ax.fill_between(N.arange(self.alpha.y.size), self.alpha.x.size * self.alpha.y, y2=self.alpha.x.size, color='white')
#		 ax.set_xlim(0, self.alpha.x.size)
#		 xticks = N.arange(N.ceil(self.alpha.x[0]), N.floor(self.alpha.x[-1]) + 1, 1) - self.alpha.x[0]
#		 xticks *= self.alpha.x.size / (self.alpha.x[-1] - self.alpha.x[0])
#		 ax.xaxis.set_ticks(xticks)
#		 def x_format(x, pos):
#			 return "%.1f" % (x * (self.alpha.x[-1] - self.alpha.x[0]) / (self.alpha.x.size) + self.alpha.x[0])
#		 ax.xaxis.set_major_formatter(FuncFormatter(x_format))
#		 yticks = N.linspace(0,1,5) * self.alpha.y.size
#		 ax.yaxis.set_ticks(yticks)
#		 def y_format(y, pos):
#			 return (y / self.alpha.y.size)
#		 ax.yaxis.set_major_formatter(FuncFormatter(y_format))
#		 ax.set_ylabel("Transmission")
#		 ax.set_xlabel("Value")
#		 pyplot.savefig(filename)

__all__ = ["ColorLinesTransferFunction"]
