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
:mod:`pymses.analysis.lic.image` --- LIC image generator module
---------------------------------------------------------------
"""
from matplotlib import cm
import numpy as N
from . import kernels
from ._line_integral_convolution import texture_streamline_convol


def lic_image(uvect, vvect, noise_texture=None, conv_kernel=None, lic_range=None, alpha=0.9, normalize=True):
    """
    Perform Line Integral Convolution upon a given vector field.

    Parameters
    ----------
    uvect: ``numpy.ndarray``
        Image x-axis vector field component 2D array.
    vvect: ``numpy.ndarray``
        Image y-axis vector field component 2D array.
    noise_texture: ``numpy.ndarray``
        noise texture 2D map of with a shape identical to the vector fields
    conv_kernel: ``numpy.ndarray``
        odd-sized 1D convolution kernel. Default None: revert to sine kernel of size 31.
    lic_range: ``list``
        lic_image value range, in [0., 1.]. Default [0.5, 0.9].
    alpha: ``float``
        Default 0.9;
    normalize: ``bool``
        Default True.
    
    
    Returns
    -------
    lic_image_rgba: ``numpy.ndarray``
        RGBA image as uint8 array of values in [0,255]
    """
    nx, ny = uvect.shape
    if lic_range is None:
        lrange = [0.5, 0.9]
    else:
        lrange = list(lic_range)
    if noise_texture is None:
        texture = N.random.rand(nx, ny).astype('f')
    else:
        texture = N.float32(noise_texture[...])
    if conv_kernel is None:
        ker = kernels.sinus_kernel().astype('f')
    else:
        ker = N.float32(conv_kernel[...])
    norm = int(normalize)

    lic_image_raw = texture_streamline_convol(uvect.astype('f'), vvect.astype('f'), texture, ker, norm)
    lrange_width = lrange[1] - lrange[0]
    lic_image_raw = (N.clip(lic_image_raw, lrange[0], lrange[1]) - lrange[0]) / lrange_width
    sm = cm.ScalarMappable(norm=None, cmap='binary')
    sm.set_clim(vmin=lrange[0], vmax=lrange[1])
    lic_data_rgba = sm.to_rgba(lic_image_raw.T, bytes=True)
    lic_data_rgba[..., 3] = N.uint8(alpha * lic_image_raw.T * 255)

    return lic_data_rgba


__all__ = ["lic_image"]
