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
    pipe.configparams.dir.ramses_output = "/data/outputs/Mass500/output_00940/"

    # Working directory (created automatically with subdirectories)
    pipe.configparams.dir.pipeline_output = "/data/postprocessing/M500_SFE01/"

The working directory will contain three subdirectories: ``ramses/`` (Zarr
archive), ``radmc3d/`` (RADMC-3D files), and ``polaris/`` (POLARIS files).


Inspecting the simulation
--------------------------

Inspect the RAMSES output before running the pipeline. Reading the file
descriptors is **not optional** for dusty simulations: it sets
``configparams.nb_dust``, the number of dust ratios / fluids that
``load_ramses`` will read. Without it the dust fields are skipped and the dust
density falls back to ``sim.dtogas`` × gas density.

.. code-block:: python

    # Hydro variables; sets nb_dust from the dust ratios
    print(pipe.read_hydro_descriptor())

    # Multi-fluid outputs only: sets nb_dust from the number of fluids
    print(pipe.read_other_file_descriptor())

    # List sinks (protostars)
    print(pipe.read_sink_info())


Configuring the AMR source
--------------------------

Select which physical fields to extract from the RAMSES output. The gas
density is always read; enable the dust field that matches your simulation:

.. code-block:: python

    pipe.configparams.amrsource.rho = True          # gas density (always on)
    pipe.configparams.amrsource.dustratios = True   # multi-grain dust ratios
    # pipe.configparams.amrsource.fluiddens = True  # multi-fluid dust densities
    pipe.configparams.amrsource.vel = False         # velocity (optional)
    pipe.configparams.amrsource.p = False           # pressure (optional)
    pipe.configparams.amrsource.temp = False        # temperature (requires pressure)

For a simulation without dust fields, leave both off: the dust density is
``sim.dtogas`` (default 0.01) times the gas density.


Defining the dust model
-----------------------

The dust model is set once in ``pipe.configparams.dust`` and used by both
POLARIS steps (opacities in Step 4, images in Step 8), so the two cannot
disagree. ``dust.mixtures`` is a dictionary of mixtures, each made of one or
more materials:

.. code-block:: python

    pipe.configparams.dust.mixtures = {
        0: {                                       # mixture id (integer)
            "silicate": {                          # free label
                "path": "/data/dust/astrosil.nk",  # optical data (.nk or tabulated)
                "distribution": "plaw",            # power law n(a) ~ a^index
                "fraction": 0.625,                 # mass fraction in the mixture
                "density": 3500.0,                 # grain bulk density [kg/m^3]
                "amin": 5e-9,                      # min grain radius [m]
                "amax": 2.5e-7,                    # max grain radius [m]
                "index": -3.5,                     # MRN slope
            },
            "graphite": {
                "path": "/data/dust/graphite.nk",
                "distribution": "plaw",
                "fraction": 0.375,
                "density": 2250.0,
                "amin": 5e-9,
                "amax": 2.5e-7,
                "index": -3.5,
            },
        },
    }

All values are SI (metres, kg/m\ :sup:`3`) and are not converted. One mixture
is needed per dust density field; for multi-grain or multi-fluid outputs a
single mixture can instead be replicated over the size bins. Every key, the
supported size distributions, and the difference between ``fraction`` and the
dust-to-gas ratio are described in :ref:`sec-dust-mixtures`.

With ``.nk`` (refractive index) files, POLARIS computes Mie cross sections at
``dust.mie_nr_sizes`` grain sizes per component (default 100):

.. code-block:: python

    pipe.configparams.dust.mie_nr_sizes = 50

Defining the observation
------------------------

Distance and wavelengths have no default and are required by Step 8:

.. code-block:: python

    pipe.configparams.imaging.distance_pc = 140.0
    pipe.configparams.imaging.wavelengths_mm = [0.87, 1.3, 3.0]
    pipe.configparams.imaging.npix = 512

``pipe.configparams.validate()`` checks the required fields before a long run.


Running the full pipeline
=========================

Step 0: Set the pymsesrc file
-----------------------------

You have the option to set, or reset the pymsesrc file according to the current RAMSES output. 
It is recommanded to set it each time you want to run a new pipeline: 

.. code-block:: python

    pipe.set_pymsesrc()

Step 1: Load RAMSES
-------------------

Read the RAMSES AMR output via pymses, flatten the octree into 1D arrays, and
store them in a compressed Zarr archive:

.. code-block:: python

    pipe.load_ramses()

This is the slowest step (dominated by pymses I/O). On subsequent runs with
the same ``pipeline_output``, you can skip it entirely: radprocess will find
the Zarr file on disk and load it directly. The same holds for every step:
see :ref:`sec-resume`.


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
for each dust mixture (i.e. each dust density field):

.. code-block:: python

    pipe.run_polaris_opacity()

The dust model is read from ``configparams.dust``. The output is a set of
``dust_mixture_XXX.dat`` files in ``polaris/data/``. Every argument can be
overridden for one call (e.g. ``dust_mixtures=...``, ``nr_threads=16``); see
:ref:`sec-params-polaris-convert`. Grain sizes are no longer set here but in
each component's ``amin``/``amax``.


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

Run POLARIS dust continuum imaging on the merged grid, with the dust model
from ``configparams.dust`` and the observation from ``configparams.imaging``:

.. code-block:: python

    pipe.render_images()

By default, images are rendered for three viewing angles (xy, xz, yz). Any
observation parameter can be overridden for one call, e.g. a single view:

.. code-block:: python

    pipe.render_images(views=["xy"], label="whole")

For zoomed-in imaging around the centre, set a field of view (full width) and a
different label so the first images are kept:

.. code-block:: python

    pipe.render_images(fov_au=2000, label="inner")   # 2000 au wide

Images are written to ``{pipeline_output}/polaris/images/{label}/{view}/``.


.. _sec-resume:

Resuming a run
==============

You never need to start from scratch. Every step writes its results to
``pipeline_output`` and the next step reads them from disk, so a run can be
stopped after any step and resumed later, in a new Python session, from the
step where it stopped. A step that fails can be fixed and re-run on its own.

What is **not** saved is the configuration: in a new session, set
``configparams`` again (at least ``dir.ramses_output``,
``dir.pipeline_output``, the dust model and the observation), then call only
the steps you need:

.. code-block:: python

    pipe = Pipeline()
    cfg = pipe.configparams
    cfg.dir.ramses_output = "/data/outputs/Mass500/output_00940/"
    cfg.dir.pipeline_output = "/data/postprocessing/M500_SFE01/"   # same as before
    cfg.dust.mixtures = {...}
    cfg.imaging.distance_pc = 140.0
    cfg.imaging.wavelengths_mm = [1.3]

    # Steps 0-3 already done: their files are in pipeline_output
    pipe.run_polaris_opacity()         # Step 4
    pipe.prepare_radmc3d_inputs()      # Step 5
    ...

Each step depends only on the outputs of earlier steps. When you change a
parameter, re-run from the first step that uses it through Step 8; the outputs
of the steps before it stay valid:

.. list-table::
   :header-rows: 1
   :widths: 46 18 36

   * - What you change
     - Re-run from
     - Files reused
   * - AMR fields (``amrsource``), ``sim.dtogas``, ``sim.use_multi_grain``,
       or a new RAMSES output
     - Step 1 (``load_ramses``)
     - none
   * - Hole around sinks (``sim.size_hole_au`` / ``hole_au``)
     - Step 2
     - Zarr store (``ramses/``)
   * - Dust model (``dust.mixtures``, ``dust.bin_edges``,
       ``dust.mie_nr_sizes``)
     - Step 4 (``run_polaris_opacity``)
     - Zarr store, both grids
   * - RADMC-3D settings (``radmc3d.*``: photons, wavelength grid, threads)
     - Step 5 (``prepare_radmc3d_inputs``)
     - Zarr store, grids, opacities
   * - Observation only (``imaging.*``: wavelengths, distance, views, field
       of view, ``npix``, ``label``)
     - Step 8 (``render_images``)
     - everything; images go to a new ``label`` folder if you change it

.. note::

   Re-running Step 4 deletes the previous POLARIS ``data/`` folder and
   ``grid_temp.dat``, so Steps 5 to 8 must follow: the temperatures and images
   would otherwise belong to the old dust model.


Complete example
================

Putting it all together:

.. code-block:: python

    from radprocess.pipeline.Pipeline import Pipeline

    pipe = Pipeline()
    pipe.configparams.dir.ramses_output = "/data/outputs/Mass500/output_00940/"
    pipe.configparams.dir.pipeline_output = "/data/postprocessing/M500_SFE01/"

    cfg = pipe.configparams

    # AMR fields to extract (sets nb_dust first)
    pipe.read_hydro_descriptor()
    cfg.amrsource.dustratios = True

    # Simulation parameters
    cfg.sim.size_hole_au = 4.0
    cfg.sim.facc = 0.5

    # Dust model (SI units). With several dust bins, this single mixture is
    # replicated over them and amin/amax come from the bins instead.
    cfg.dust.mixtures = {
        0: {
            "silicate": {"path": "/data/dust/astrosil.nk", "distribution": "plaw",
                         "fraction": 0.625, "density": 3500.0,
                         "amin": 5e-9, "amax": 2.5e-7, "index": -3.5},
            "graphite": {"path": "/data/dust/graphite.nk", "distribution": "plaw",
                         "fraction": 0.375, "density": 2250.0,
                         "amin": 5e-9, "amax": 2.5e-7, "index": -3.5},
        },
    }

    # Observation
    cfg.imaging.distance_pc = 140.0
    cfg.imaging.wavelengths_mm = [0.87, 1.3, 3.0]
    cfg.imaging.npix = 512

    cfg.validate()

    # Run the pipeline
    pipe.set_pymsesrc()                # Step 0
    pipe.load_ramses()                 # Step 1
    pipe.convert_to_polaris()          # Step 2
    pipe.convert_to_radmc()            # Step 3
    pipe.run_polaris_opacity()         # Step 4
    pipe.prepare_radmc3d_inputs()      # Step 5
    pipe.run_radmc3d_mctherm()         # Step 6
    pipe.merge_temperature()           # Step 7
    pipe.render_images()               # Step 8


Configuration reference
=======================

All pipeline parameters are accessible through ``pipe.configparams``, in the
groups ``dir``, ``amrsource``, ``sim``, ``dust``, ``radmc3d``, ``polaris`` and
``imaging``. Each field, its type, default, unit and the method arguments that
override it are listed in :ref:`chap-parameters`. The configuration is strict:
assigning a field that does not exist raises an error immediately.

To inspect all current values:

.. code-block:: python

    print(pipe.configparams)
    print(pipe.configparams.radmc3d)
    print(pipe.configparams.dust)