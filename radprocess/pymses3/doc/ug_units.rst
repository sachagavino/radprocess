Dealing with units
==================

.. topic:: Need

    Okay, I have read my data quite easily. What are the units of these data ? How do I convert
    them into human-readable units ?
    
    **Example** : From a RAMSES hydro simulation, I want to convert my density field unit into
    the H/cc unit.



Dimensional physical constants
------------------------------

In pymses, a specific module has been designed for this purpose : :mod:`~pymses.utils.constants`.

It contains a bunch of useful dimensional physical constants (expressed in ISU) which you can use for unit conversion factors computation, adimensionality tests, etc. ::

    from pymses.utils import constants as C
    print C.kpc
    Out:(3.085677e+19 m)
    print C.Msun
    Out:(1.9889e+30 kg)

Each constant is an :class:`~pymses.utils.constants.Unit` instance, on which you can call the :meth:`~pymses.utils.constants.Unit.express` method to convert this constant into
another dimension-compatible constant. If the dimensions are not compatible, a ValueError will be raised ::

    factor = C.kpc.express(C.ly)
    print "1 kpc = %f ly"%factor
    Out:1 kpc = 3261.563163 ly

    print C.Msun.express(C.km)
    Out:ValueError: Incompatible dimensions between (1.9889e+30 kg) and (1000.0 m)


Basic operations between these constants are enabled ::

    unit_density = 1.0E9 * C.Msun / C.kpc**3
    print "1Msun/kpc**3 = %f H/cc"%unit_density.express(C.H_cc)
    Out:1Msun/kpc**3 = 30.993246 H/cc


RAMSES data units
-----------------

The units of each RAMSES output data are read from the output info file. You can manipulate the values of these units by using the *info* parameter (see :ref:`ramses_output`) ::

    ro = RamsesOutput("/data/simu/outputs", 10)

    ro.info
    Out:{'H0': 1.0,
        'aexp': 1.0,
        'boxlen': 200.0,
        'dom_decomp': <pymses.sources.ramses.hilbert.HilbertDomainDecomp object at 0x9df0aac>,
        'levelmax': 14,
        'levelmin': 7,
        'ncpu': 64,
        'ndim': 3,
        'ngridmax': 1000000,
        'nstep_coarse': 1200,
        'omega_b': 0.0,
        'omega_k': 0.0,
        'omega_l': 0.0,
        'omega_m': 1.0,
        'time': 10.724764558171801,
        'unit_density': (6.77025430199e-20 m^-3.kg),
        'unit_length': (6.17135516256e+21 m),
        'unit_mass': (1.9891e+39 kg),
        'unit_pressure': (2.91283226304e-10 m^-1.kg.s^-2),
        'unit_temperature': (517290.92492 K),
        'unit_time': (4.70430312424e+14 s),
        'unit_velocity': (65592.6605874 m.s^-1)}


Assuming you already have sampled the AMR density field of this output into a *pdset* :class:`~pymses.core.datasets.PointDataset` containing all your sampling points (see :ref:`amr_sampling`), you can convert your density field (in code unit) into the unit of your choice::

    rho_field_H_cc = pdset["rho"] * ro.info["unit_density"].express(C.H_cc)





.. topic:: Warning

    You must take care of manipulating RAMSES data in an unit-coherent way !!!

        - **Good**::

            info = ro.info

            # Density
            rho_H_cc = rho_ramses * info["unit_density"].express(C.H_cc)

            # Mass
            part_mass_Msun = part_mass * info["unit_mass"].express(C.Msun)
 
            # Kinetic energy
            factor = (info["unit_mass"] * info["unit_velocity"]**2).express(C.J)
            kin_energy_J = part_mass * part_vel**2 * factor
        
        - **Not so good**::
        
            info = ro.info

            # Density
            factor = (info["unit_mass"] / info["unit_length"]**3).express(C.H_cc)
            rho_H_cc = rho_ramses * factor

            # Mass
            factor = (info["unit_density"]*info["unit_length"]**3).express(C.Msun)
            part_mass_Msun = part_mass * factor
 
            # Kinetic energy
            factor = (info["unit_pressure"] * info["unit_length"]**3).express(C.J)
            kin_energy_J = part_mass * part_vel**2 * factor

