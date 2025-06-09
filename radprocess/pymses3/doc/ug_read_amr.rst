AMR data access
###############

AMR data source
***************

If you want to deal with the AMR data, you need to call the :meth:`~pymses.sources.ramses.output.RamsesOutput.amr_source` method of your :class:`~pymses.sources.ramses.output.RamsesOutput` object with a single argument which is a list of the different fields you might need in your analysis.

When calling the :meth:`~pymses.sources.ramses.output.RamsesOutput.amr_source`, the fields you have access to are :

  * "rho" : the gas density field
  * "vel" : the gas velocity field
  * "P" : the gas pressurre field
  * "g" : the gravitational acceleration field

To modify the list of available data fields, see :doc:`ug_formats`.

.. code-block:: python

    ro = pymses.RamsesOutput("/data/Aquarius/output", 193)
    amr = ro.amr_source(["rho", "vel", "P", "g"])


.. topic:: Warning
    
    The data source you just created does not contain data. It is designed to provide an *on-demand* access to the data. To be memory-friendly, nothing is read from the disk yet at this point. All the :file:`amr_00193.out_*`, :file:`hydro_00193.out_*` and :file:`grav_00193.out_*` files are only linked to the data source for further processing.


AMR data handling
*****************

AMR data is a bit more complicated to handle than particle data. To perform various analysis, PyMSES provides you with two different tools to get your AMR data :

  * :ref:`cells_to_points`
  * :ref:`amr_sampling`


.. _cells_to_points:


AMR grid to cell list conversion
================================

The :class:`~pymses.filters.CellsToPoints` filter converts the AMR tree structure into a :class:`~pymses.core.datasets.IsotropicExtPointDataset` containing a list of the AMR grid `leaf envelope` cells :

    * The *points* parameter of the datasets coming from the generated data source will contain the coordinates of the cell centers.
    * These datasets will have an additional :meth:`~pymses.core.datasets.IsotropicExtPointDataset.get_sizes` method giving you the size of each cell.

.. code-block:: python
    
    from pymses.filters import CellsToPoints
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    [...]
    # Cell centers
    ccenters = cells.points
    # Cell sizes
    dx = cells.get_sizes()

.. topic:: Warning

    As a :class:`Filter<pymses.core.sources.Filter>`,  the `cell_source` object you first created is
    another data provider, it doesn't contain actual data. To read the data, use
    :meth:`~pymses.sources.ramses.sources.RamsesGenericSource.get_domain_dset`,
    :meth:`~pymses.core.sources.DataSource.iter_dsets` or
    :meth:`~pymses.core.sources.DataSource.flatten` method as described in :doc:`ug_read_particles`.




.. _amr_sampling:

AMR field point-sampling
========================

Another way to read the AMR data is to perform a sampling of the AMR fields with a set of sampling points coordinates of your choice.
In PyMSES, this is done quite easily  with the :func:`~pymses.analysis.point_sampler.sample_points` function :

.. code-block:: python

    from pymses.analysis import sample_points
    sample_dset = sample_points(amr, points)

  
The returned `sample_dset` will be a :class:`~pymses.core.datasets.PointDataset` containing all your sampling
points and the corresponding value of the different AMR fields you selected.


.. topic:: Note

    In backstage, the point sampling is performed with a `tree search` algorithm, which makes this particular process of AMR data access both **user-friendly** and **efficient**.


For example, this method can be used :

 * for visualization purposes (see :doc:`ug_visu_slices`).
 * when computing profiles (see :doc:`ug_ana_profiles`)

