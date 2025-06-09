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
:mod:`pymses.analysis.lic.kernels` --- LIC kernels module
---------------------------------------------------------

Implementation of :
 - Basic sinus kernel function
 - Hanning ripples according to B. Cabral and L. Leedom (1993) to animate LIC images.
"""
import numpy as N


def sinus_kernel(kernel_length=31):
    ker = N.sin(N.arange(kernel_length)*N.pi/kernel_length)
    return ker / N.sum(ker)


def __hanning_window(x):
    """
    Return the value of the hanning function.

    One cycle corresponds to the domain [0,1].
     
    Parameters
    ----------
    x: ``list``
        x value list.
    """
    x = N.array(x)
    return 0.5 - 0.5 * N.cos(2.0 * N.pi * x)


def hanning_ripples(kernel_length=31, shift=0, ripples=3):
    """
    Return a kernel composed of an hanning envelope superposed on phase shifted hanning ripples. 

    A vector field can be "animated" to provide a sense of the direction of flow by creating a series of lic images,
    each one computed from the same texture, but with a sequence of kernels that are slightly phase shifted with
    respect to each other. 

    Parameters
    ----------
    kernel_length : int
      Kernel length. 
    shift : float
      Phase shift [0,1].
    ripples : int
      Number of ripples. 

    Returns
    -------
    out : 1D array (N,)
      Linearly spaced kernel values. 
    """

    """
    Example (CHECK THAT IT WORKS)
    -------
    vectors = ...
    texture = ...
    images = []
    L = 10
    for i in range(L):
       kernel = hanning_ripples(shift=N.linspace(0,1,L,endpoint=False)).astype(N.float32)
       images.append(line_integral_convolution(vectors, texture, kernel)
    """

    enveloppe = N.hanning(kernel_length)
    signal = __hanning_window(shift + N.linspace(0, 1, kernel_length) * ripples)
    return signal * enveloppe


__all__ = ['sinus_kernel', 'hanning_ripples']
