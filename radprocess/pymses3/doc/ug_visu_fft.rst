FFT-convolved maps
===================

Intro
-----

A very simple, fast and accurate data projection (3D->2D) method : each particle/AMR cell is convolved by a 2D gaussian kernel (`Splatter`) which size depends on the local AMR grid level.

The convolution of the binned particles/AMR cells histogram with the gaussian kernels is performed with FFT techniques by a :class:`~pymses.analysis.visualization.fft_projection.MapFFTProcessor`. You can see two examples of this method below :

 * :ref:`particles_fft_map`
 * :ref:`amr_fft_map`

.. Topic:: Important note on operators

    You must keep in mind that any `X` :class:`~pymses.analysis.operator.AbstractOperator` you use with this method must describe an **extensive** physical variable since this method compute a summation over particle/AMR quantities :

    :math:`map[i,j] = \displaystyle\sum_{\text{particles/AMR cells}} X`

Examples
--------

.. _particles_fft_map:

Particles map
*************

.. plot:: pyplots/fft_map_part.py
    :include-source:

.. _amr_fft_map:

AMR data map
************

.. plot:: pyplots/fft_map_amr.py
    :include-source:
