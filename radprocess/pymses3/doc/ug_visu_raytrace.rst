Ray-traced maps
===============

Intro
-----

Ray-traced maps are computed in PyMSES by integrating a physical quantity along `rays`, each one corresponding to a pixel of the map. Ray-tracing is handled by a :class:`~pymses.analysis.visualization.raytracing.RayTracer`. You can see two examples of this method below :

 * :ref:`ray_trace_density`
 * :ref:`ray_trace_Temp_min`
 * :ref:`ray_trace_levelmax`

.. Topic:: Important note on operators

    You must keep in mind that any `X` :class:`~pymses.analysis.operator.AbstractOperator` you use with this method must describe an **intensive** physical variable since this method compute an integral of an AMR quantity over each pixel surface and along the line-of-sight :

    :math:`map[i,j] = \displaystyle\int_{z_{\text{min}}}^{z_{\text{max}}} X \textrm{d}S_{\text{pix}}\textrm{d}z`

Examples
--------

.. _ray_trace_density:

Density map
***********

.. plot:: pyplots/ray_trace_rho.py
    :include-source:


.. _ray_trace_Temp_min:

Min. temperature map
********************

.. plot:: pyplots/ray_trace_Tmin.py
    :include-source:


.. _ray_trace_levelmax:

Max. AMR level of refinement map
********************************

.. plot:: pyplots/ray_trace_lmax.py
    :include-source:

Multiprocessing
---------------

    If you are using python 2.6 or higher, the RayTracer will try to use multiprocessing speed up. You can desactivate it to save RAM memory and processor use by setting the multiprocessing option to False:

.. code-block:: python

    map = rt.process(scal_op, cam, multiprocessing = False)