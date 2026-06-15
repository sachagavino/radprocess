.. _chap-subboxes:

Working with subboxes
*********************

When a RAMSES simulation contains multiple sinks (protostars), the full-domain
pipeline processes the entire cloud at once. For studying individual
protostellar disks, you would prefer to use **subboxes**: small octree grids extracted around each sink
and processed independently. This saves much more of CPU-time than performing the radiative transfer on the full grid,
in particular if one only needs to generate images centered on a single sink.

This page walks through the complete subbox workflow, from extraction to
science analysis.


Overview
========

The subbox pipeline differs from the full-domain pipeline in several ways:

- **Extraction**: one subbox per sink, centered on the sink position, with a
  padded grid for octree alignment.
- **Temperature**: RADMC-3D ``mctherm`` runs independently per sink (parallelisable).
- **Imaging**: POLARIS renders each subbox with a field of view and detector
  shift centered on the sink.
- **Analysis**: spectral indices, optically-thin masses, and bolometric
  luminosities are computed per sink.

The grid is extracted larger than the requested field of view (padded) to
ensure artifact-free octree alignment. The imaging tools then crop to the
user's requested FOV centered on the sink. This means the star in
``stars.inp`` is at an offset from the grid center, that is by design.


Setting up
==========

Start with the same pipeline setup as the basic workflow:

.. code-block:: python

    from radprocess.pipeline.Pipeline import Pipeline

    pipe = Pipeline()
    pipe.configparams.dir.ramses_output = "/data/outputs/Mass500/output_00940/"
    pipe.configparams.dir.pipeline_output = "/data/postprocessing/M500_subboxes/"

    # Define the AMR fields you want to extract
    pipe.configparams.amrsource.rho = True
    pipe.configparams.amrsource.dust = True

    # Simulation parameters
    pipe.configparams.sim.size_hole_au = 4.0
    pipe.configparams.sim.facc = 0.1

    # Dust materials
    dust = [
        {"path": "/data/dust/silicate_d03.nk", "weight": 0.625},
        {"path": "/data/dust/graphite.nk",     "weight": 0.375},
    ]

    # Load RAMSES (Step 1, same as basic pipeline)
    pipe.load_ramses()


Step 2: Extract subboxes
========================

Extract subboxes around each sink in both POLARIS and RADMC-3D format:

.. code-block:: python

    # POLARIS format
    pipe.convert_subboxes(
        which_rad="polaris",
        box_half_width_au=500,       # ±500 AU around each sink
        isolation_radius_au=200,     # skip sinks closer than 200 AU to a neighbour
        min_cells=1000,              # select boxes with a minimum of cell numbers. 
        hole_au=0.0,                 # dig hole around sink (this replaces the value set in the Simulation parameters)
    )

    # RADMC-3D format (same extraction, different output format)
    pipe.convert_subboxes(
        which_rad="radmc",
        box_half_width_au=500,
    )

This creates one folder per sink under ``polaris/subboxes/`` and
``radmc3d/subboxes/``. Each folder contains the grid files and metadata:

- ``sink_offset.txt``: the sink's offset from the grid center (in metres for
  POLARIS, cm for RADMC-3D).
- ``requested_hw_au.txt``: the user's requested half-width (for automatic FOV
  in imaging).

.. note::

    The actual subbox is **larger** than the requested ±500 AU. The code pads
    the box to ensure grid-aligned octree boundaries, then crops to your
    requested FOV during imaging. The log output shows the actual box size and
    the margin to the box edge.


Step 3: Run POLARIS opacity
===========================

Same as the basic pipeline — this only needs to run once:

.. code-block:: python

    pipe.run_polaris_opacity(dust_components=dust)


Step 4: Prepare RADMC-3D inputs
===============================

Prepare wavelength grid, opacities, and per-sink ``stars.inp`` and
``radmc3d.inp`` files:

.. code-block:: python

    pipe.prepare_radmc3d_inputs(
        subbox=True,
        nphot=5_000_000,
        n_wavelengths=200,
        wave_min=0.27,       # µm (match the dustkappa range)
        wave_max=3000,       # µm
        scattering_mode=2,
    )

This distributes the shared files (opacities, wavelength grid) to all subbox
folders and writes a per-sink ``stars.inp`` (star at the sink offset) and
``radmc3d.inp`` (including ``subbox_*`` parameters for the regrid, centered on
the sink).


Step 5: Run RADMC-3D mctherm
=============================

Run the thermal Monte Carlo computation for each sink. This can be done
sequentially:

.. code-block:: python

    pipe.run_radmc3d_mctherm(subbox="all")

Or for a single sink:

.. code-block:: python

    pipe.run_radmc3d_mctherm(subbox="sink_0026")


Running in parallel with HTCondor
---------------------------------

For large numbers of sinks, run ``mctherm`` in parallel using HTCondor.
Create three files:

**run_mctherm_single.py**:

.. code-block:: python

    import sys
    from radprocess.pipeline.Pipeline import Pipeline

    pipe_path = "/data/postprocessing/M500_subboxes/"
    radmc3d_path = "/path/to/radmc3d"

    sink_name = sys.argv[1]

    pipe = Pipeline()
    pipe.configparams.dir.pipeline_output = pipe_path
    pipe.run_radmc3d_mctherm(subbox=sink_name, radmc3d_binary=radmc3d_path)

**jobs_mctherm.sh**:

.. code-block:: bash

    #!/bin/bash
    set -e
    source /path/to/radprocess_env/bin/activate
    SINK_NAME=$1
    python /path/to/run_mctherm_single.py ${SINK_NAME}

**jobs_mctherm.sub**:

.. code-block:: text

    universe   = vanilla
    executable = /path/to/jobs_mctherm.sh
    arguments  = $(sink_name)
    initialdir = /data/postprocessing/M500_subboxes/
    output     = std_outputs/mctherm-$(sink_name).out
    error      = std_outputs/mctherm-$(sink_name).err
    log        = std_outputs/mctherm-$(sink_name).log
    should_transfer_files = NO
    request_memory = 31GB
    request_cpus   = 14
    queue sink_name from sink_list.txt

Generate the sink list and submit:

.. code-block:: bash

    ls -d radmc3d/subboxes/sink_* | xargs -n1 basename > sink_list.txt
    chmod +x jobs_mctherm.sh
    condor_submit jobs_mctherm.sub

This launches one job per sink, all running in parallel.


Step 6: Merge temperature
=========================

Inject the RADMC-3D temperatures into the POLARIS grids:

.. code-block:: python

    pipe.merge_temperature(subbox="all")


Step 7: Render images
=====================

Render POLARIS dust emission images for each sink:

.. code-block:: python

    pipe.render_images(
        dust_components=dust,
        npix=256,
        distance_pc=140.0,
        wavelengths_mm=[0.87, 1.3, 3.0],
        subbox="all",
        fov_au=1000,     # 1000 AU wide image (±500 AU, matching the extraction)
        polaris_binary="/path/to/polaris",
    )

The ``fov_au`` parameter sets the field of view in AU. It should match the
``box_au`` used for the RADMC-3D regrid so that density and intensity maps are
directly comparable. The image is automatically centered on the sink via the
POLARIS detector shift.

You can also render POLARIS images for a single or a set of chosen sinks:

.. code-block:: python

    pipe.render_images(
        dust_components=dust,
        npix=256,
        distance_pc=140.0,
        wavelengths_mm=[0.87, 1.3, 3.0],
        subbox=["sink_0026", "sink_0040"],
        fov_au=1000,     # 1000 AU wide image (±500 AU, matching the extraction)
        polaris_binary="/path/to/polaris",
    )

For computing bolometric luminosities, render with a broader wavelength range:

.. code-block:: python

    pipe.render_images(
        dust_components=dust,
        npix=128,
        distance_pc=140.0,
        wavelengths_mm=[
            0.001, 0.002, 0.005, 0.01, 0.02, 0.05,   # 1–50 µm
            0.07, 0.1, 0.16, 0.25, 0.45,               # 70–450 µm
            0.87, 1.3, 3.0,                             # mm
        ],
        subbox="all",
        fov_au=1000,
        polaris_binary="/path/to/polaris",
    )


Diagnostic plots
================

Density mosaic
--------------

After running ``radmc3d subbox_regrid`` in each sink folder, plot the
projected dust column density:

.. code-block:: python

    from radprocess.plotting import plot

    fig = plot.subbox_density_mosaic(
        pipe_path + "radmc3d/subboxes/",
        axis="z",
        log=True,
        cmap="inferno",
    )

Intensity mosaic
----------------

Plot the POLARIS Stokes I images:

.. code-block:: python

    fig = plot.subbox_mosaic(
        pipe_path,
        quantity="intensity",
        view="xy",
        wavelength_idx=0,      # first wavelength (0.87 mm)
        log=True,
        cmap="inferno",
    )

Other quantities from the same FITS files:

.. code-block:: python

    # Optical depth at 1.3 mm
    fig = plot.subbox_mosaic(pipe_path, quantity="tau", view="xy", wavelength_idx=1)

    # Column density
    fig = plot.subbox_mosaic(pipe_path, quantity="column_density", view="xy")


Science analysis
================

Bolometric luminosity and temperature
-------------------------------------

Compute T_bol and L_bol for each sink from the multi-wavelength SED:

.. code-block:: python

    table = plot.compute_all_tbol_lbol(
        pipe_path,
        ramses_path="/data/outputs/Mass500/output_00940/",
        view="xy",
        distance_pc=140.0,
        aperture_au=200,
    )

    # BLT diagram
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.scatter(table["T_bol"], table["L_bol_obs"])
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$T_{\rm bol}$ [K]")
    ax.set_ylabel(r"$L_{\rm bol,obs}$ [$L_\odot$]")

    # Compare observed vs true luminosity
    fig2, ax2 = plt.subplots()
    ax2.scatter(table["L_bol_true"], table["L_bol_obs"])
    ax2.plot([1e-3, 1e3], [1e-3, 1e3], 'k--', alpha=0.3)
    ax2.set_xscale("log"); ax2.set_yscale("log")
    ax2.set_xlabel(r"$L_{\rm bol,true}$ [$L_\odot$]")
    ax2.set_ylabel(r"$L_{\rm bol,obs}$ [$L_\odot$]")


Spectral index and optically-thin mass
---------------------------------------

Compute the mm spectral index α and compare optically-thin mass estimates
with the true mass within the same aperture:

.. code-block:: python

    table = plot.compute_alpha_mass(
        pipe_path,
        ramses_path="/data/outputs/Mass500/output_00940/",
        view="xy",
        distance_pc=140.0,
        wavelengths_mm=[0.87, 1.3, 3.0],
        aperture_au=200,
        T_dust_K=20.0,            # observer's assumed temperature
        kappa_ref_cm2g=2.3,       # Beckwith+1990 (includes gas-to-dust=100)
        kappa_ref_mm=1.3,
        kappa_beta=1.7,
        gas_to_dust_ratio=100.0,
    )

The returned table contains matched arrays for each sink:

- ``alpha``: spectral index from a power-law fit to the 3 mm fluxes
- ``M_thin_Msun``: optically-thin mass at each wavelength (shape ``N×3``)
- ``M_true_aperture_Msun``: true (gas+dust) mass within the aperture
- ``M_sink_Msun``: sink particle mass from RAMSES
- ``M_ratio``: M_thin / M_true — the mass recovery fraction

.. code-block:: python

    # α vs mass recovery
    fig, ax = plt.subplots()
    ax.scatter(table["alpha"], table["M_ratio"])
    ax.set_xlabel(r"$\alpha_{\rm mm}$")
    ax.set_ylabel(r"$M_{\rm thin} / M_{\rm true}$")
    ax.axhline(1, ls='--', color='k', alpha=0.3)

    # Optically-thin mass vs true mass
    fig2, ax2 = plt.subplots()
    ax2.scatter(table["M_true_aperture_Msun"], table["M_thin_Msun"][:, 1])
    ax2.plot([1e-4, 10], [1e-4, 10], 'k--', alpha=0.3)
    ax2.set_xscale("log"); ax2.set_yscale("log")
    ax2.set_xlabel(r"$M_{\rm true}$ (aperture) [$M_\odot$]")
    ax2.set_ylabel(r"$M_{\rm thin}$ (1.3 mm) [$M_\odot$]")

Points below the 1:1 line indicate optically thick sources where the mass is
underestimated. Lower α correlates with higher optical depth and worse mass
recovery.