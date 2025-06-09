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


class ConvolKernel(object):
    """Convolution kernel class
    """

    def __init__(self, ker_func, size_func=None, max_size=None):
        """Convolution kernel builder

        ker_func  : convolution kernel function => 2D function lambda x, y, size: f(x, y, size)
        size_func : kernel size factor. The size of the convolution kernel is set to 'size_func' x the local cell size.
        max_size  : maximum size of the convolution kernel.
        """
        # Kernel 2D function
        self.ker_func = ker_func

        # Size function
        if size_func is None:
            self.size_func = lambda dset: dset.get_sizes()
        else:
            self.size_func = size_func

        # Max. size of the convolution kernel
        self.max_size = max_size

        # Allow to change the size factor without having to create a new kernel object
        self.FFTkernelSizeFactor = 1

    def get_size(self, dset):
        """
        """
        return self.FFTkernelSizeFactor * self.size_func(dset)

    def get_max_size(self):
        return self.FFTkernelSizeFactor * self.max_size

    def get_convolved_axes(self):
        raise NotImplementedError

    def convol_fft(self, kernel_size, map_dict, ext_camera, verbose=None):
        """
        FFT convolution method designed to convolve a dict. of maps into a single map.

        Parameters
        ----------
        kernel_size: ``float``
            size of the comvolution kernel
        map_dict : ``dict``
            map dict. where the dict. keys are the size of the convolution kernel.
        ext_camera : ExtendedCamera
            Camera corrsponding to maps of the map dict.
        verbose: ``bool``
            verbosity boolean flag. Default None.

        Returns
        -------
        TODO
        """
        map = dict.fromkeys(list(map_dict.keys()))

        # Map coordinates edges
        xc, yc = ext_camera.get_pixels_coordinates_centers()

        # Region of interest limits
        imin, imax, jmin, jmax = ext_camera.get_window_mask()

        level = int(N.log2(1.0 / kernel_size))
        if verbose is None or verbose:
            print("   * level = %i" % level)
        kernel = self.ker_func(xc, yc, kernel_size)
        if N.max(kernel) != 0.0:
            kernel = N.fft.fftshift(kernel / N.sum(kernel))
            kernel = N.fft.fftn(kernel)

            for var in map_dict:
                conv = kernel * N.fft.fftn(map_dict[var])
                dat = N.real(N.fft.ifftn(conv))
                del conv
                map[var] = dat[imin:imax, jmin:jmax]
                del dat
        else:
            if verbose is None or verbose:
                print("WARNING : this kernel is too small to be taken into account on this map.")

        return map


class GaussSplatter2DKernel(ConvolKernel):
    """2D Gaussian splatter convolution kernel
    """

    def __init__(self, size_func=None, max_size=None):
        """2D Gaussian splatter convolution kernel builder
        """
        def ker_func(x, y, size):
            sigma = size / 2.
            ker = 1 / (2 * N.pi * sigma ** 2) * N.outer(N.exp(-x ** 2 / (2.0 * sigma ** 2)),
                                                        N.exp(-y ** 2 / (2.0 * sigma ** 2)))
            return ker

        super(GaussSplatter2DKernel, self).__init__(ker_func, size_func, max_size)

    def get_convolved_axes(self):
        return [0, 1]


class Cos2SplatterKernel(ConvolKernel):
    """2D Squared cosine splatter convolution kernel
    """

    def __init__(self, size_func=None, max_size=None):
        """2D Squared cosine splatter convolution kernel builder
        """

        def f(x, y, size):
            ker = N.outer(N.cos(N.pi * x / (2. * size)) ** 2, N.cos(N.pi * y / (2. * size)) ** 2)
            ker[(N.abs(x) > size), :] = 0.0
            ker[:, (N.abs(y) > size)] = 0.0
            return ker

        super(Cos2SplatterKernel, self).__init__(f, size_func, max_size)

    def get_convolved_axes(self):
        return [0, 1]


class PyramidSplatterKernel(ConvolKernel):
    """2D pyramidal splatter convolution kernel
    """

    def __init__(self, size_func=None, max_size=None):
        """2D pyramidal splatter convolution kernel builder
        """

        def f(x, y, size):
            ker = N.outer(N.abs(1.0 - x / size), N.abs(1.0 - y / size))
            ker[(N.abs(x) > size), :] = 0.0
            ker[:, (N.abs(y) > size)] = 0.0
            return ker

        super(PyramidSplatterKernel, self).__init__(f, size_func, max_size)

    def get_convolved_axes(self):
        return [0, 1]


class GaussSplatter1DKernel(ConvolKernel):
    """
    1D Gaussian splatter convolution kernel
    """

    def __init__(self, axis, size_func=None, max_size=None):
        assert axis in [0, 1], "axis param must be in [0, 1]."
        self.axis = axis
        if axis == 0:
            def ker_func(x, y, size):
                yvect = N.zeros_like(y)
                yvect[(N.abs(y) == N.min(N.abs(y)))] = 1.0
                sigma = size/2.
                ker = N.outer(1 / (N.sqrt(2*N.pi) * sigma) * N.exp(-x ** 2 / (2.0 * sigma ** 2)), yvect)
                return ker
        else:
            def ker_func(x, y, size):
                xvect = N.zeros_like(x)
                xvect[(N.abs(x) == N.min(N.abs(x)))] = 1.0
                sigma = size/2.
                ker = N.outer(xvect, 1 / (N.sqrt(2*N.pi) * sigma) * N.exp(-y ** 2 / (2. * sigma ** 2)))
                return ker
        super(GaussSplatter1DKernel, self).__init__(ker_func, size_func, max_size)

    def get_convolved_axes(self):
        return [self.axis]


__all__ = ["ConvolKernel", "GaussSplatter2DKernel", "GaussSplatter1DKernel", "PyramidSplatterKernel",
           "Cos2SplatterKernel"]
