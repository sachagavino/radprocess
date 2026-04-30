.. _chap-basicuse:

How to use radprocess
*********************

This section walks through a complete pipeline run, from loading a RAMSES
simulation to producing synthetic dust continuum images. Each step is shown
as it would appear in a Jupyter notebook.


Setting up the pipeline
=======================

Start by importing the ``Pipeline`` class and creating an instance:

.. code-block:: python

    from radprocess.pipeline.Pipeline import Pipeline

    pipe = Pipeline()


Pointing to the data
--------------------

Two paths must be set before anything else: the RAMSES output directory and a
working directory where radprocess will write all intermediate and final files.

.. code-block:: python

    # RAMSES simulation output
    pipe.configparams.ramsesoutput.ramses_output_dir = "/data/outputs/Mass500/output_00940/"

    # Working directory (created automatically with subdirectories)
    pipe.set_working_dir("/data/postprocessing/M500_SFE01")

The working directory will contain three subdirectories: ``ramses/`` (Zarr
archive), ``radmc3d/`` (RADMC-3D files), and ``polaris/`` (POLARIS files).


Inspecting the simulation
--------------------------

Before running the pipeline, you can inspect the RAMSES output:

.. code-block:: python

    # Check available hydro variables
    print(pipe.read_hydro_descriptor())

    # List sinks (protostars)
    print(pipe.read_sink_info())


Configuring the AMR source
--------------------------

Select which physical fields to extract from the RAMSES output. At minimum,
you need gas density and dust ratios enabled:

.. code-block:: python

    pipe.configparams.amrsource.rho = True       # gas density
    pipe.configparams.amrsource.dust = True       # dust density ratios
    pipe.configparams.amrsource.vel = False       # velocity (optional)
    pipe.configparams.amrsource.p = False         # pressure (optional)
    pipe.configparams.amrsource.temp = False      # temperature (requires pressure)


Defining dust components
------------------------

The dust material properties are needed by both POLARIS (Steps 4 and 8) and
are passed as a list of dictionaries. Each entry points to a POLARIS-format
cross-section file (``.cs``) and a mass fraction weight:

.. code-block:: python

    dust_components = [
        {"path": "/data/dust/silicate_d03.cs", "weight": 0.625},
        {"path": "/data/dust/carbon_z96.cs",   "weight": 0.375},
    ]

For a single-component model (e.g., pure silicate), use a single entry with
``weight = 1.0``.


Running the full pipeline
=========================

Step 1: Load RAMSES
-------------------

Read the RAMSES AMR output via pymses, flatten the octree into 1D arrays, and
store them in a compressed Zarr archive:

.. code-block:: python

    pipe.load_ramses()

This is the slowest step (dominated by pymses I/O). On subsequent runs, you
can skip it entirely: radprocess will find the Zarr file on disk and load it
directly.


Step 2: Convert to POLARIS grid
-------------------------------

Build a POLARIS-format binary octree from the Zarr data:

.. code-block:: python

    pipe.convert_to_polaris()

The output is a single binary file (``ramses_grid_XXXXX.dat``) in the
``polaris/`` subdirectory. Densities are converted from CGS (g/cm\ :sup:`3`)
to SI (kg/m\ :sup:`3`) as required by POLARIS.


Step 3: Convert to RADMC-3D grid
---------------------------------

Build the RADMC-3D octree grid and dust density files:

.. code-block:: python

    pipe.convert_to_radmc()

This writes ``amr_grid.inp``, ``dust_density.inp``, and ``stars.inp`` into
the ``radmc3d/`` subdirectory, all in CGS units.


Step 4: Run POLARIS opacity
----------------------------

Execute POLARIS with a single photon package to generate dust opacity tables
for each size bin:

.. code-block:: python

    pipe.run_polaris_opacity(dust_components=dust_components)

The output is a set of ``dust_mixture_XXX.dat`` files in ``polaris/data/``.
You can override the default grain size range and power-law index either by
passing them directly or by setting them in the configuration:

.. code-block:: python

    # Via the config dataclass (persistent for the session)
    pipe.configparams.polaris.dust_size_min = 5e-9       # 5 nm
    pipe.configparams.polaris.dust_size_max = 2.5e-7     # 250 nm
    pipe.configparams.polaris.dust_size_powerlaw = -3.5   # MRN

    # Or via keyword arguments (one-off override)
    pipe.run_polaris_opacity(
        dust_components=dust_components,
        dust_size_min=1e-8,
        dust_size_max=1e-5,
    )


Step 5: Prepare RADMC-3D inputs
--------------------------------

Convert the POLARIS opacity tables to RADMC-3D format and write all remaining
input files:

.. code-block:: python

    pipe.prepare_radmc3d_inputs()

This writes ``dustkappa_*.inp``, ``dustopac.inp``, ``wavelength_micron.inp``,
and ``radmc3d.inp`` into the ``radmc3d/`` subdirectory. The control parameters
for ``radmc3d.inp`` can be adjusted through the configuration:

.. code-block:: python

    pipe.configparams.radmc3d.nphot = 10_000_000
    pipe.configparams.radmc3d.setthreads = 16
    pipe.configparams.radmc3d.scattering_mode = 1
    pipe.configparams.radmc3d.modified_random_walk = 1


Step 6: Run RADMC-3D mctherm
-----------------------------

Execute the RADMC-3D Monte Carlo thermal computation:

.. code-block:: python

    temp_file = pipe.run_radmc3d_mctherm()

This runs ``radmc3d mctherm`` in the ``radmc3d/`` directory and produces
``dust_temperature.bdat``. The computation time depends on the number of
photon packages and the optical depth of the model.


Step 7: Merge temperature
--------------------------

Inject the RADMC-3D dust temperatures back into the POLARIS grid:

.. code-block:: python

    pipe.merge_temperature()

This produces ``grid_temp.radmc3d.dat`` in the ``polaris/`` subdirectory. The
original POLARIS grid is preserved as ``grid_temp.polaris.dat``.


Step 8: Render images
----------------------

Run POLARIS dust continuum imaging on the merged grid:

.. code-block:: python

    pipe.render_images(
        dust_components=dust_components,
        npix=512,
        distance_pc=140.0,
        wavelengths_mm=[0.87, 1.3, 3.0],
    )

By default, images are rendered for three viewing angles (xy, xz, yz). You can
select specific views:

.. code-block:: python

    pipe.render_images(
        dust_components=dust_components,
        npix=512,
        distance_pc=140.0,
        wavelengths_mm=[0.87, 1.3, 3.0],
        views=["xy"],
        label="whole",
    )

For zoomed-in imaging around the central region, set a field of view and
increase the midplane zoom:

.. code-block:: python

    pipe.render_images(
        dust_components=dust_components,
        npix=512,
        distance_pc=140.0,
        wavelengths_mm=[0.87, 1.3, 3.0],
        midplane_zoom=10,
        fov_m=1e14,        # field of view half-width in metres
        label="inner",
    )

Images are written to ``{working_dir}/images/{label}/{view}/``.


Complete example
================

Putting it all together:

.. code-block:: python

    from radprocess.pipeline.Pipeline import Pipeline

    pipe = Pipeline()
    pipe.configparams.ramsesoutput.ramses_output_dir = "/data/outputs/Mass500/output_00940/"
    pipe.set_working_dir("/data/postprocessing/M500_SFE01")

    # AMR fields to extract
    pipe.configparams.amrsource.rho = True
    pipe.configparams.amrsource.dust = True

    # Simulation parameters
    pipe.configparams.sim.size_hole_au = 4.0
    pipe.configparams.sim.facc = 0.5

    # Dust materials
    dust = [
        {"path": "/data/dust/silicate_d03.cs", "weight": 0.625},
        {"path": "/data/dust/carbon_z96.cs",   "weight": 0.375},
    ]

    # Run the pipeline
    pipe.load_ramses()                                      # Step 1
    pipe.convert_to_polaris()                               # Step 2
    pipe.convert_to_radmc()                                 # Step 3
    pipe.run_polaris_opacity(dust_components=dust)           # Step 4
    pipe.prepare_radmc3d_inputs()                            # Step 5
    pipe.run_radmc3d_mctherm()                               # Step 6
    pipe.merge_temperature()                                 # Step 7
    pipe.render_images(                                      # Step 8
        dust_components=dust,
        npix=512,
        distance_pc=140.0,
        wavelengths_mm=[0.87, 1.3, 3.0],
    )


Configuration reference
=======================

All pipeline parameters are accessible through ``pipe.configparams``. The main
groups are:

``pipe.configparams.sim``
    Simulation-level parameters: ``size_hole_au`` (hole radius around sinks),
    ``facc`` (accretion efficiency factor).

``pipe.configparams.radmc3d``
    RADMC-3D control parameters: ``nphot``, ``nphot_scat``, ``setthreads``,
    ``scattering_mode``, ``modified_random_walk``, ``rto_style``,
    ``rto_single``, ``wave_min``, ``wave_max``, ``n_wavelengths``.

``pipe.configparams.polaris``
    POLARIS parameters: ``dust_size_min``, ``dust_size_max``,
    ``dust_size_powerlaw``, ``mean_molecular_weight``, ``mass_fraction``,
    ``nr_threads``, ``polaris_binary``.

``pipe.configparams.amrsource``
    Flags for which AMR fields to extract: ``rho``, ``dust``, ``vel``, ``p``,
    ``temp``, ``bl``, ``br``.

To inspect all current values:

.. code-block:: python

    print(pipe.configparams)
    print(pipe.configparams.radmc3d)
    print(pipe.configparams.polaris)