Get a RAMSES output into PyMSES
###############################

.. topic:: Use case

    You want to select a specific RAMSES output directory and get somme basic information about it


.. _ramses_output:

RAMSES output selection
***********************

First, you need to select the snapshot of your RAMSES simulation you would like to read by creating a :class:`~pymses.sources.ramses.output.RamsesOutput` object :

.. code-block:: python

    import pymses
    ro = pymses.RamsesOutput("/data/Aquarius/outputs", 193)

In this example, you are intersted in the files contained in :file:`/data/Aquarius/output/ouput_00193/`

.. _info:

Ouput information
*****************

To get some details about this specific output/simulation. Everything you need is in the `info` parameter :

.. code-block:: python

    ro.info
    Out:{'H0': 73.0,
        'aexp': 1.0000620502295701,
        'boxlen': 1.0,
        'dom_decomp': <pymses.sources.ramses.hilbert.HilbertDomainDecomp object at 0x3305e10>,
        'levelmax': 18,
        'levelmin': 7,
        'ncpu': 64,
        'ndim': 3,
        'ngridmax': 800000,
        'nstep_coarse': 9578,
        'omega_b': 0.039999999105930301,
        'omega_k': 0.0,
        'omega_l': 0.75,
        'omega_m': 0.25,
        'time': 6.2446534480863097e-05,
        'unit_density': (2.50416381926e-27 m^-3.kg),
        'unit_length': (4.21943976727e+24 m),
        'unit_mass': (1.88116596007e+47 kg),
        'unit_pressure': (2.50385294276e-13 m^-1.kg.s^-2),
        'unit_temperature': (12021826243.9 K),
        'unit_time': (4.21970170037e+17 s),
        'unit_velocity': (9999379.26156 m.s^-1)}
    
    ro.info["ncpu"]
    Out:64

    ro.info["boxlen"] / 2**ro.info["levelmax"]
    Out:3.814697265625e-06

This way, you can easily find the units of your data (see :doc:`ug_units`).


