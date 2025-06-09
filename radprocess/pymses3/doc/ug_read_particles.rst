Reading particles
#################

Particle data source
********************

If you want to look at the particles, you need to create a :class:`~pymses.sources.ramses.sources.RamsesParticleSource`.
To do so, call the :meth:`~pymses.sources.ramses.output.RamsesOutput.particle_source` method of your :class:`~pymses.sources.ramses.output.RamsesOutput` object with a list of the different fields you might need in your analysis.

The available fields are :

  * "vel" : the velocities of the particles
  * "mass" : the mass of the particles
  * "id" : the id number of the particles
  * "level" : the AMR level of refinement of the cell each particle belongs to
  * "epoch" : the birth time of the particles (0.0 for ic particles, >0.0 for star formation particles)
  * "metal" : the metallicities of the particles


.. code-block:: python

    ro = pymses.RamsesOutput("/data/Aquarius/output", 193)
    part = ro.particle_source(["vel", "mass"])



.. topic:: Warning
    
    The data source you just created does not contain data. It is designed to provide an *on-demand* access to the data. To be memory-friendly, nothing is read from the disk yet at this point. All the :file:`part_00193.out_*` files are only linked to the data source for further processing.


.. _point_dataset:

PointDataset
************

At the opposite, a :class:`~pymses.core.datasets.PointDataset` is an actual data container. 


Single CPU particle dataset
===========================

If you want to read all the particles of the cpu number 3 (written in :file:`part_00193.out_00003`), use the :meth:`~pymses.sources.ramses.sources.RamsesGenericSource.get_domain_dset` method :


.. code-block:: python

    dset3 = part.get_domain_dset(3)
    Out:Reading particles : /data/Aquarius/output/output_00193/part_00193.out00003


Number of particles
===================

Every :class:`~pymses.core.datasets.PointDataset` has a *npoints* ``int`` parameter which gives you the number of particles in this dataset :

.. code-block:: python

    print "CPU 3 has %i particles"%dset3.npoints
    Out:CPU 3 has 157976 particles

Particle coordinates
====================

The *points* parameter of the :class:`~pymses.core.datasets.PointDataset` contains the coordinates of the particles :

.. code-block:: python

    print dset3.points
    Out:array([[ 0.49422911,  0.51383241,  0.50130034],
            [ 0.49423128,  0.51374527,  0.50136899],
            [ 0.49420231,  0.51378629,  0.50190981],
            ..., 
            [ 0.49447162,  0.51394969,  0.50146777],
            [ 0.49422794,  0.51378071,  0.50176276],
            [ 0.4946566 ,  0.51491008,  0.50117673]])


Particle fields
===============

You also have an easy access to the different fields :

.. code-block:: python

    print dset3["mass"]
    Out:array([  4.69471978e-07,   4.69471978e-07,   9.38943957e-07, ...,
                4.69471978e-07,   4.69471978e-07,   4.69471978e-07])



Whole point data source concatenation
*************************************

To read all the particles from all the ncpus :file:`part_00193.out*` files and concatenate them into a single (but maybe not memory-friendly) dataset, call the :meth:`~pymses.core.sources.DataSource.flatten` method of your *part* object :


.. code-block:: python

    dset_all = part.flatten()
    Out:Reading particles : /data/Aquarius/output/output_00193/part_00193.out00001
        Reading particles : /data/Aquarius/output/output_00193/part_00193.out00002
        Reading particles : /data/Aquarius/output/output_00193/part_00193.out00003
        Reading particles : /data/Aquarius/output/output_00193/part_00193.out00004
        
        [...]
    
        Reading particles : /data/Aquarius/output/output_00193/part_00193.out00062
        Reading particles : /data/Aquarius/output/output_00193/part_00193.out00063
        Reading particles : /data/Aquarius/output/output_00193/part_00193.out00064

    print "Domain has %i particles"%dset_all.npoints
    Out:Domain has 10000000 particles


CPU-by-CPU particles
********************

In most cases, you won't have enough memory to load all the particles of your simulation domain into a single dataset. You have two different options :

  * Filter your particles (see :doc:`ug_data_filtering`).
  * Your analysis can be done on a cpu-by-cpu basis. The :class:`~pymses.sources.ramses.sources.RamsesParticleSource` provides a :meth:`~pymses.core.sources.DataSource.iter_dsets` iterator yielding cpu-by-cpu datasets :

.. code-block:: python

    for dset in part.iter_dsets():
        print dset.npoints

    Out:Reading particles : /data/Aquarius/output/output_00193/part_00193.out00001
        254210
        Reading particles : /data/Aquarius/output/output_00193/part_00193.out00002
        214330
        Reading particles : /data/Aquarius/output/output_00193/part_00193.out00003
        359648
        [...]
        Reading particles : /data/Aquarius/output/output_00193/part_00193.out00064
        351203


