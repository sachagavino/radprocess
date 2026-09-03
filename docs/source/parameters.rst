.. _chap-parameters:

Parameters reference
********************************************

This page gathers every parameter a user can set when driving the Radprocess
pipeline, grouped into four stages:

#. :ref:`sec-params-ramses` — reading and extracting the RAMSES output.
#. :ref:`sec-params-polaris-convert` — converting the RAMSES output to POLARIS.
#. :ref:`sec-params-radmc-convert` — converting the RAMSES output to RADMC-3D.
#. :ref:`sec-params-render` — rendering synthetic images with POLARIS.

.. note::

   Parameters reach the pipeline through **two** mechanisms:

   * **Configuration parameters** — persistent fields stored on the
     configuration object (``pipe.configparams``). They are set once and
     reused by every stage. These are the fields defined in
     :mod:`radprocess.utils.config`.
   * **Method arguments** — passed at call time to a pipeline method
     (``convert_to_polaris``, ``render_images``, ...). Many method arguments
     default to ``None``; when left as ``None`` they fall back to the matching
     configuration parameter.

Setting configuration parameters
================================

.. code-block:: python

   from radprocess.pipeline import Pipeline

   pipe = Pipeline()
   cfg = pipe.configparams          # a ConfigParams instance

   # Stage 1 — RAMSES extraction
   cfg.dir.ramses_output = "ramses_outputs/"
   cfg.sim.size_hole_au  = 5.0
   cfg.amrsource.vel     = True

   # Stage 2 / 4 — POLARIS
   cfg.polaris.dust_size_max = 1.0e-3

   # Stage 3 — RADMC-3D
   cfg.radmc3d.nphot = 2_000_000

The configuration object is *strict*: assigning an unknown field, or a value of
the wrong type, raises an error. It also renders as a rich table in a Jupyter
notebook and as a formatted table in the terminal.


.. _sec-params-ramses:

1. Conversion from RAMSES
=========================

Controls how the RAMSES output is located, which fields are read, and how the
domain is extracted (central hole, sub-boxes) before any radiative-transfer
conversion. These settings are shared by both the POLARIS and RADMC-3D paths.

Directories (configuration)
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 22 10 24 44

   * - Parameter
     - Type
     - Default
     - Description
   * - ``ramses_output``
     - str
     - ``'ramses_outputs/'``
     - The RAMSES simulation output directory path.
   * - ``pipeline_output``
     - str
     - ``'pipeline_outputs/'``
     - The main pipeline output directory path.

Source fields (configuration)
-----------------------------

Which AMR fields are read from the RAMSES output. ``rho`` is always included;
the remaining fields are read only if enabled *and* present in the output.

.. list-table::
   :header-rows: 1
   :widths: 18 8 10 64

   * - Parameter
     - Type
     - Default
     - Description
   * - ``rho``
     - bool
     - ``True``
     - Include the gas density field (always True).
   * - ``fluiddens``
     - bool
     - ``False``
     - Include the multiple fluid densities (if they exist).
   * - ``dustratios``
     - bool
     - ``False``
     - Include the multiple dust ratios (if they exist).
   * - ``vel``
     - bool
     - ``False``
     - Include the velocity field.
   * - ``fluidvel``
     - bool
     - ``False``
     - Include the fluid velocity field (if they exist).
   * - ``bl``
     - bool
     - ``False``
     - Include the left magnetic field (B) components.
   * - ``br``
     - bool
     - ``False``
     - Include the right magnetic field (B) components.
   * - ``p``
     - bool
     - ``False``
     - Include the pressure field.
   * - ``xi``
     - bool
     - ``False``
     - Include the ionization fraction field.
   * - ``phi``
     - bool
     - ``False``
     - Include the gravitational potential field.
   * - ``g``
     - bool
     - ``False``
     - Include the gravitational acceleration field.
   * - ``temp``
     - bool
     - ``False``
     - Include the gas temperature field.

Simulation / physics (configuration)
------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 24 8 12 56

   * - Parameter
     - Type
     - Default
     - Description
   * - ``size_hole_au``
     - float
     - ``0.0``
     - [AU] Size of the central hole to extrude around the sinks.
   * - ``dtogas``
     - float
     - ``0.01``
     - Dust-to-gas mass ratio. Used only if this is not multi-grain mode, or if no dust ratios exist in the RAMSES output.
   * - ``facc``
     - float
     - ``0.1``
     - Accretion fraction converted into radiation. Change at your own risk!
   * - ``use_ramses_T``
     - bool
     - ``True``
     - Use RAMSES stellar temperatures (if available in the sink info) as stellar inputs for the RT simulations.
   * - ``use_ramses_acc_rate``
     - bool
     - ``True``
     - Use RAMSES accretion rates (if available in the sink info) to derive stellar radii for the RT simulations.
   * - ``use_multi_grain``
     - bool
     - ``True``
     - Do RT with multiple bins if True (when the RAMSES output has multiple bins). If False, dust density is computed as ``dtogas * rho``.
   * - ``nb_dust``
     - int
     - ``0``
     - Number of dust species / bins (bookkeeping; populated during conversion).

Extraction & grid controls (method arguments)
---------------------------------------------

Passed to ``convert_to_radmc``, ``convert_to_polaris`` and
``convert_subboxes``. Grid-layout options (``gridstyle``, ``coordsystem``)
apply to RADMC-3D output only.

.. list-table::
   :header-rows: 1
   :widths: 24 16 14 46

   * - Argument
     - Type
     - Default
     - Description
   * - ``gridstyle``
     - str
     - ``"octtree"``
     - Grid layout for RADMC-3D output: ``"octtree"`` or ``"regular"``.
   * - ``coordsystem``
     - str
     - ``"cartesian"``
     - Coordinate system for RADMC-3D output.
   * - ``hole_au``
     - float or None
     - ``None``
     - Hole radius around sinks (AU). If None, falls back to ``sim.size_hole_au``.
   * - ``box_half_width_au``
     - float
     - ``100.0``
     - Half-width of each sub-box in AU (100 gives 200 AU boxes).
   * - ``isolation_radius_au``
     - float
     - ``100.0``
     - Minimum separation for sink filtering, in AU.
   * - ``boxlen_pc``
     - float or None
     - ``None``
     - RAMSES box size in pc. If None, derived from the Zarr store.
   * - ``require_luminosity``
     - bool
     - ``True``
     - If True, skip sinks without luminosity data.
   * - ``sink_indices``
     - list of int or None
     - ``None``
     - If provided, skip filtering and process only these sinks.
   * - ``min_cells``
     - int
     - ``0``
     - Minimum number of cells in a sub-box; smaller sinks are skipped (0 = no minimum).
   * - ``which_rad``
     - str
     - ``"radmc"``
     - Sub-box output format: ``"radmc"`` or ``"polaris"``.


.. _sec-params-polaris-convert:

2. RAMSES output converted to POLARIS
=====================================

Settings used when writing the POLARIS octree grid and generating the dust
opacity tables (``convert_to_polaris`` and ``run_polaris_opacity``). This stage
also inherits the extraction controls from Stage 1.

POLARIS dust & runtime (configuration)
--------------------------------------

These fields live in ``cfg.polaris`` and are shared with the rendering stage
(Stage 4).

.. list-table::
   :header-rows: 1
   :widths: 26 8 14 52

   * - Parameter
     - Type
     - Default
     - Description
   * - ``dust_size_min``
     - float
     - ``5e-9``
     - [m] Minimum grain radius.
   * - ``dust_size_max``
     - float
     - ``2.5e-7``
     - [m] Maximum grain radius.
   * - ``dust_size_powerlaw``
     - float
     - ``-3.5``
     - Power-law exponent for the grain size distribution (e.g. -3.5 for MRN).
   * - ``mean_molecular_weight``
     - float
     - ``2.37``
     - Mean molecular weight (mu) of the gas.
   * - ``mass_fraction``
     - float
     - ``1.0``
     - Dust-to-gas mass fraction.
   * - ``nr_threads``
     - int
     - ``8``
     - Number of OpenMP threads for POLARIS.
   * - ``polaris_binary``
     - str
     - ``'polaris'``
     - Name or path of the POLARIS executable.

Opacity run (method arguments)
------------------------------

Passed to ``run_polaris_opacity``. Arguments left as ``None`` fall back to the
``cfg.polaris`` fields above.

.. list-table::
   :header-rows: 1
   :widths: 26 20 12 42

   * - Argument
     - Type
     - Default
     - Description
   * - ``dust_mixtures``
     - dict
     - *(required)*
     - Dust material definitions (paths, size distribution, fractions). Must be supplied explicitly since it contains paths specific to the user's setup.
   * - ``mean_molecular_weight``
     - float or None
     - ``None``
     - Gas mean molecular weight. Falls back to ``polaris.mean_molecular_weight``.
   * - ``mass_fraction``
     - float or None
     - ``None``
     - Dust-to-gas mass fraction. Falls back to ``polaris.mass_fraction``.
   * - ``nr_threads``
     - int or None
     - ``None``
     - OpenMP threads. Falls back to ``polaris.nr_threads``.
   * - ``grid_path``
     - str/Path or None
     - ``None``
     - Path to the POLARIS grid file. If None, auto-detected.
   * - ``n_dust_override``
     - int or None
     - ``None``
     - Override the auto-detected number of dust species.
   * - ``polaris_binary``
     - str or None
     - ``None``
     - POLARIS executable. Falls back to ``polaris.polaris_binary``.
   * - ``cleanup``
     - bool
     - ``True``
     - If True, remove previous POLARIS run outputs before starting.

The grid-writing step ``convert_to_polaris`` additionally accepts ``hole_au``
(see Stage 1).


.. _sec-params-radmc-convert:

3. RAMSES output converted to RADMC-3D
======================================

Settings used when converting POLARIS opacities to RADMC-3D format, writing the
RADMC-3D input files, and running the thermal Monte-Carlo
(``prepare_radmc3d_inputs``, ``thermal_radmc3d``, ``run_radmc3d_mctherm``).
This stage also inherits the grid controls from Stage 1
(``convert_to_radmc``: ``gridstyle``, ``coordsystem``).

RADMC-3D settings (configuration)
---------------------------------

Fields in ``cfg.radmc3d``.

.. list-table::
   :header-rows: 1
   :widths: 26 8 14 52

   * - Parameter
     - Type
     - Default
     - Description
   * - ``nphot``
     - int
     - ``1000000``
     - Number of photon packages for mctherm.
   * - ``nphot_scat``
     - int
     - ``1000000``
     - Number of photon packages for scattering Monte Carlo.
   * - ``setthreads``
     - int
     - ``8``
     - Number of OpenMP threads for RADMC-3D.
   * - ``scattering_mode``
     - int
     - ``1``
     - Scattering mode (1 = isotropic, 2 = anisotropic with HG, 5 = full).
   * - ``scattering_mode_max``
     - int
     - ``1``
     - Maximum scattering mode.
   * - ``modified_random_walk``
     - int
     - ``1``
     - Enable modified random walk (1 = yes, 0 = no). Accelerates optically thick regions.
   * - ``rto_style``
     - int
     - ``3``
     - Output style for dust_temperature (1 = ascii, 3 = binary).
   * - ``rto_single``
     - int
     - ``1``
     - Single-precision output (1 = yes, 0 = no).
   * - ``wave_min``
     - float
     - ``0.27``
     - [micron] Minimum wavelength for the wavelength grid.
   * - ``wave_max``
     - float
     - ``3000.0``
     - [micron] Maximum wavelength for the wavelength grid.
   * - ``n_wavelengths``
     - int
     - ``200``
     - Number of log-spaced wavelength points.

Input preparation (method arguments)
-------------------------------------

Passed to ``prepare_radmc3d_inputs``. Arguments left as ``None`` fall back to
the ``cfg.radmc3d`` fields above.

.. list-table::
   :header-rows: 1
   :widths: 24 18 12 46

   * - Argument
     - Type
     - Default
     - Description
   * - ``polaris_data_dir``
     - str/Path or None
     - ``None``
     - Path to the POLARIS ``data/`` directory. If None, auto-detected.
   * - ``n_dust``
     - int or None
     - ``None``
     - Number of dust species. If None, auto-detected from the POLARIS output.
   * - ``wave_min``
     - float or None
     - ``None``
     - [micron] Wavelength grid minimum. Falls back to ``radmc3d.wave_min``.
   * - ``wave_max``
     - float or None
     - ``None``
     - [micron] Wavelength grid maximum. Falls back to ``radmc3d.wave_max``.
   * - ``n_wavelengths``
     - int or None
     - ``None``
     - Number of wavelength points. Falls back to ``radmc3d.n_wavelengths``.
   * - ``nphot``
     - int or None
     - ``None``
     - Photon packages for mctherm. Falls back to ``radmc3d.nphot``.
   * - ``setthreads``
     - int or None
     - ``None``
     - OpenMP threads. Falls back to ``radmc3d.setthreads``.
   * - ``scattering_mode``
     - int or None
     - ``None``
     - Scattering mode. Falls back to ``radmc3d.scattering_mode``.
   * - ``subbox``
     - bool
     - ``False``
     - If True, distribute the shared input files to all sub-box folders.

Thermal Monte-Carlo (method arguments)
--------------------------------------

``thermal_radmc3d`` writes the RADMC-3D input files and (optionally) runs the
thermal Monte-Carlo. The ``write_*`` flags toggle individual input files.

.. list-table::
   :header-rows: 1
   :widths: 24 10 12 54

   * - Argument
     - Type
     - Default
     - Description
   * - ``run``
     - bool
     - ``True``
     - Run the thermal Monte-Carlo. If False, assume the RADMC-3D output files already exist.
   * - ``nphot``
     - int/float
     - ``1e4``
     - Number of photon packages for this thermal run.
   * - ``write_opac``
     - bool
     - ``True``
     - Write the dust opacity input files.
   * - ``write_control``
     - bool
     - ``True``
     - Write the ``radmc3d.inp`` control file.
   * - ``write_star``
     - bool
     - ``True``
     - Write the stellar source file.
   * - ``write_wave``
     - bool
     - ``True``
     - Write the wavelength grid file.
   * - ``write_mcmono``
     - bool
     - ``True``
     - Write the monochromatic Monte-Carlo wavelength file.
   * - ``write_ext``
     - bool
     - ``True``
     - Write the external / additional input files.

Running mctherm (method arguments)
----------------------------------

Passed to ``run_radmc3d_mctherm``.

.. list-table::
   :header-rows: 1
   :widths: 24 22 14 40

   * - Argument
     - Type
     - Default
     - Description
   * - ``radmc3d_binary``
     - str
     - ``"radmc3d"``
     - Name or path of the RADMC-3D executable.
   * - ``subbox``
     - None / str / list / "all"
     - ``None``
     - Which directory to run in: None (full cloud), a single sink, a list of sinks, or "all" sub-boxes.


.. _sec-params-render:

4. Rendering images with POLARIS
================================

Settings used to merge the RADMC-3D temperature into the POLARIS grid and to
produce the synthetic dust-continuum images
(``merge_temperature`` and ``render_images``). Most parameters here are method
arguments; a few default to ``cfg.polaris`` fields (Stage 2).

Merging the temperature (method arguments)
------------------------------------------

Passed to ``merge_temperature``.

.. list-table::
   :header-rows: 1
   :widths: 20 24 14 42

   * - Argument
     - Type
     - Default
     - Description
   * - ``n_dust``
     - int or None
     - ``None``
     - Number of dust species. If None, auto-detected.
   * - ``subbox``
     - None / str / list / "all"
     - ``None``
     - Which grid to merge the temperature into.

Rendering (method arguments)
----------------------------

Passed to ``render_images``.

.. list-table::
   :header-rows: 1
   :widths: 24 22 16 38

   * - Argument
     - Type
     - Default
     - Description
   * - ``dust_mixtures``
     - dict
     - *(required)*
     - Dust material definitions (same structure as the opacity run).
   * - ``npix``
     - int
     - *(required)*
     - Image resolution (npix x npix pixels).
   * - ``distance_pc``
     - float
     - *(required)*
     - Source distance in parsecs.
   * - ``wavelengths_mm``
     - list of float
     - *(required)*
     - Wavelengths to image, in millimetres.
   * - ``views``
     - list of str or None
     - ``None``
     - Viewing angles to render. None renders all three standard views: ``xy``, ``xz``, ``yz``.
   * - ``fov_m``
     - float or None
     - ``None``
     - Field of view in metres (full width). Overridden by ``fov_au`` if set. None = full grid extent.
   * - ``fov_au``
     - float or None
     - ``None``
     - Field of view in AU (full width). For sub-boxes, use the same value as the RADMC-3D regrid box.
   * - ``label``
     - str
     - ``"whole"``
     - Output subdirectory label (e.g. "whole" or "inner").
   * - ``grid_path``
     - str/Path or None
     - ``None``
     - Path to the merged grid. If None, auto-detected.
   * - ``n_dust``
     - int or None
     - ``None``
     - Number of dust species. If None, auto-detected.
   * - ``nr_threads``
     - int or None
     - ``None``
     - OpenMP threads. Falls back to ``polaris.nr_threads``.
   * - ``mean_molecular_weight``
     - float or None
     - ``None``
     - Gas mu. Falls back to ``polaris.mean_molecular_weight``.
   * - ``mass_fraction``
     - float or None
     - ``None``
     - Dust-to-gas mass fraction. Falls back to ``polaris.mass_fraction``.
   * - ``polaris_binary``
     - str or None
     - ``None``
     - POLARIS executable. Falls back to ``polaris.polaris_binary``.
   * - ``cleanup_views``
     - bool
     - ``True``
     - Remove previous image output for each view before rendering.
   * - ``polaris_cmd``
     - str
     - ``"CMD_DUST_EMISSION"``
     - POLARIS command: ``CMD_DUST_EMISSION`` (thermal), ``CMD_DUST_SCATTERING`` (scattered light), or both in sequence.
   * - ``alignment``
     - str
     - ``"ALIG_PA"``
     - Grain alignment mechanism: ``ALIG_PA``, ``ALIG_IDG``, ``ALIG_RAT``, ``ALIG_INTERNAL``, or empty/None for none.
   * - ``peel_off``
     - bool
     - ``True``
     - Use the peel-off technique for scattering (more efficient for images; only relevant for scattering).
   * - ``acceptance_angle``
     - float or None
     - ``None``
     - Acceptance angle for scattered light, in degrees. None uses the POLARIS default.
   * - ``nr_photons_scat``
     - int or None
     - ``None``
     - Number of photon packages for the scattering Monte-Carlo.
   * - ``source_star_scat``
     - list of dict or None
     - ``None``
     - Stellar sources for scattering (each with position, radius, temperature).
   * - ``subbox``
     - None / str / list / "all"
     - ``None``
     - Which target to render: None (full cloud), a single sink, a list, or "all" sub-boxes.
