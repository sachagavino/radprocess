Data filtering
==============

In PyMSES, a ``Filter`` is a data source that :
  * filter the data coming from an origin data source.
  * provides a new data source out of this filtering process.


Region filter
-------------

For a lot of analysis, you are often interested in a particular region of your simulation domain, for example :

  - spherical region centered on a dark matter halo in a cosmological simulation.
  - cylindrical region containing a galactik disk or a cosmological filament.
  - ...

.. code-block:: python

    # Region of interest
    from pymses.utils.regions import Sphere
    center = [0.5, 0.5, 0.5]
    radius = 0.1
    region = Sphere(center, radius)


To filter data source with a region, use the :class:`~pymses.filters.RegionFilter`::

    from pymses.filters import RegionFilter    
    from pymses import RamsesOutput
    ro = RamsesOutput("/data/Aquarius/output/", 193)
    parts = ro.particle_source(["mass"])
    amr = ro.amr_source(["rho"])

    # Particle filtering
    filt_parts = RegionFilter(region, parts)

    # AMR data filtering
    filt_amr = RegionFilter(region, amr)


.. topic:: Note

    The region filters can greatly improve the I/O performance of your analysis process since it doesn't require to read all the cpu files (of your entire simulation domain) but only those whose domain intersects your region of interest.

    The filtering process occurs not only at the cpu level (as any other :class:`~pymses.core.sources.Filter`) but also in the choice of required cpu files.


The CellsToPoints filter
------------------------

see :ref:`cells_to_points`.


Function filters
----------------

You can also define your own function filter. Here an example where only the particles of mass equal to :math:`3\times10^{3} M_{\odot}` are gathered :

.. code-block:: python

    from pymses.filters import PointFunctionFilter 
    from pymses.utils import constants as C
    
    part_source = ro.particle_source(["mass"])

    # Stellar disc particles filter : only keep particles of mass = 3000.0 Msun 
    part_mass_Msun = 3.0E3 * C.Msun
    part_mass_code = part_mass_Msun.express(ro.info["unit_mass"])
    st_disc_func = lambda dset: (dset["mass"]==part_mass_code)

    # Stellar disc particle data source
    st_disc_parts = PointFunctionFilter(st_disc_func, part_source)


Randomly decimated data
-----------------------

You can use the :class:`~pymses.filters.PointRandomDecimatedFilter` filter to pick up only a given fraction of points (randomly chosen) from your point-based data::

    from pymses.filters import PointRandomDecimatedFilter
    part_source = ro.particle_source(["mass"])

    # Pick up 10 % of the particles
    fraction = 0.1
    dec_parts = PointRandomDecimatedFilter(fraction, part_source)


Combining filters
-----------------

You can pile up as many filters as you want to get the particular data you're interested in::

    # Region filter
    reg_parts = RegionFilter(region, parts)

    # Stellar disc filter
    st_disc_parts = PointFunctionFilter(st_disc_func, reg_parts)

    # 10 % randomly decimated filter
    dec_parts = PointRandomDecimatedFilter(fraction, st_disc_parts)

In this example, the `dec_parts` data source will provide you 10% of the stellar particles contained within a given `region`
