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
     :mod:`radprocess.utils.config`, grouped into ``dir``, ``amrsource``,
     ``sim``, ``dust``, ``radmc3d``, ``polaris`` and ``imaging``.
   * **Method arguments** — passed at call time to a pipeline method
     (``convert_to_polaris``, ``render_images``, ...). Config-backed method
     arguments default to ``None``; when left as ``None`` they are resolved
     from the matching configuration parameter, so an explicit call value
     always wins and an omitted one inherits ``configparams``.

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

   # Stage 2 / 4 — dust setup (shared by opacity and imaging)
   cfg.dust.mixtures          = {...}     # see DustConfig
   cfg.dust.mean_molecular_weight = 2.37
   cfg.dust.mass_fraction     = 0.01

   # Stage 3 — RADMC-3D
   cfg.radmc3d.nphot = 2_000_000

   # Stage 4 — imaging
   cfg.imaging.distance_pc    = 140.0
   cfg.imaging.wavelengths_mm = [1.3]

Before a run, ``cfg.validate()`` checks the required fields
(``dust.mixtures``, ``imaging.distance_pc`` and ``imaging.wavelengths_mm``).
The configuration object is *strict*: assigning an unknown field, or a value
of the wrong type, raises an error. It also renders as a rich table in a
Jupyter notebook and as a formatted table in the terminal.


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
Enabling ``temp`` implicitly forces ``p`` on, since temperature is derived from
pressure.

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
   * - ``use_multi_grain``
     - bool
     - ``True``
     - Do RT with multiple bins if True (when the RAMSES output has multiple bins). If False, dust density is computed as ``dtogas * rho``, ignoring any per-bin data.
   * - ``use_ramses_T``
     - bool
     - ``True``
     - *Reserved — not currently wired.* Intended to use RAMSES stellar temperatures as stellar inputs. See the note below.
   * - ``use_ramses_acc_rate``
     - bool
     - ``True``
     - *Reserved — not currently wired.* Intended to use RAMSES accretion data to derive stellar radii. See the note below.

.. warning::

   ``use_ramses_T`` and ``use_ramses_acc_rate`` are currently **inactive**:
   they can be set but are not read anywhere, so changing them has no effect.
   The stellar temperature is always taken from the RAMSES ``Teff`` column and
   the radius is derived from the total (internal + accretion) luminosity. The
   intended semantics — in particular whether the accretion *rate* (with
   ``facc``) or the accretion *luminosity* should drive the radius — are still
   under discussion, so these flags are left unwired on purpose.

.. note::

   ``nb_dust`` is a top-level ``ConfigParams`` field (not part of ``sim``); it
   is bookkeeping for the number of dust species and is populated during
   conversion.

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

Dust composition (configuration)
--------------------------------

The dust setup lives in ``cfg.dust`` and is the single source of truth shared
by both the opacity step (Stage 2) and rendering (Stage 4).

.. list-table::
   :header-rows: 1
   :widths: 26 10 14 50

   * - Parameter
     - Type
     - Default
     - Description
   * - ``mixtures``
     - dict
     - ``{}``
     - Dust mixtures, ``{mixture_id: {component: {path, distribution, fraction, density, amin, amax, index}}}``. Required (grain sizes come from each component's ``amin``/``amax``/``distribution``).
   * - ``mean_molecular_weight``
     - float
     - ``2.37``
     - Mean molecular weight (mu) of the gas.
   * - ``mass_fraction``
     - float
     - ``1``
     - Gas-to-dust mass ratio: 100 already applied during grid conversion (written as the POLARIS ``<mass_fraction>`` keyword).

.. note::

   Within a mixture, each component's ``fraction`` is its mass fraction *within
   the dust* (the components of a mixture sum to 1). This is distinct from
   ``mass_fraction`` above, which is the overall dust-to-gas ratio.

POLARIS execution (configuration)
---------------------------------

.. list-table::
   :header-rows: 1
   :widths: 26 8 14 52

   * - Parameter
     - Type
     - Default
     - Description
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
``cfg.dust`` and ``cfg.polaris`` fields above.

.. list-table::
   :header-rows: 1
   :widths: 26 20 12 42

   * - Argument
     - Type
     - Default
     - Description
   * - ``dust_mixtures``
     - dict or None
     - ``None``
     - Dust material definitions. If None, uses ``dust.mixtures``.
   * - ``mean_molecular_weight``
     - float or None
     - ``None``
     - Gas mean molecular weight. Falls back to ``dust.mean_molecular_weight``.
   * - ``mass_fraction``
     - float or None
     - ``None``
     - Dust-to-gas mass fraction. Falls back to ``dust.mass_fraction``.
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
(``prepare_radmc3d_inputs`` and ``run_radmc3d_mctherm``). This stage also
inherits the grid controls from Stage 1 (``convert_to_radmc``: ``gridstyle``,
``coordsystem``).

RADMC-3D settings (configuration)
---------------------------------

Fields in ``cfg.radmc3d``. All of these are now forwarded through
``prepare_radmc3d_inputs`` to the single ``radmc3d.inp`` writer.

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
     - Number of photon packages for the scattering Monte-Carlo.
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
   * - ``nphot_scat``
     - int or None
     - ``None``
     - Scattering photon packages. Falls back to ``radmc3d.nphot_scat``.
   * - ``setthreads``
     - int or None
     - ``None``
     - OpenMP threads. Falls back to ``radmc3d.setthreads``.
   * - ``scattering_mode``
     - int or None
     - ``None``
     - Scattering mode. Falls back to ``radmc3d.scattering_mode``.
   * - ``scattering_mode_max``
     - int or None
     - ``None``
     - Maximum scattering mode. Falls back to ``radmc3d.scattering_mode_max``.
   * - ``modified_random_walk``
     - int or None
     - ``None``
     - Modified random walk. Falls back to ``radmc3d.modified_random_walk``.
   * - ``rto_style``
     - int or None
     - ``None``
     - Output style. Falls back to ``radmc3d.rto_style``.
   * - ``rto_single``
     - int or None
     - ``None``
     - Single-precision output. Falls back to ``radmc3d.rto_single``.
   * - ``subbox``
     - bool
     - ``False``
     - If True, distribute the shared input files to all sub-box folders.

.. note::

   The ``radmc3d.inp`` control file is written by a single writer,
   ``radmc3d.prepare.write_radmc3d_control``. Beyond the explicit keywords
   above it accepts an ``extra_options`` dict for any other ``radmc3d.inp``
   key (used, for example, to pass the ``subbox_*`` regrid cube per sink).

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


.. _sec-params-render:

4. Rendering images with POLARIS
================================

Settings used to merge the RADMC-3D temperature into the POLARIS grid and to
produce the synthetic dust-continuum images (``merge_temperature`` and
``render_images``). The observation is defined by ``cfg.imaging``; the dust
setup is inherited from ``cfg.dust`` (Stage 2) and execution settings from
``cfg.polaris``.

Imaging / observation (configuration)
-------------------------------------

Fields in ``cfg.imaging``. ``distance_pc`` and ``wavelengths_mm`` have no safe
default and are required (checked by ``validate()``).

.. list-table::
   :header-rows: 1
   :widths: 28 12 16 44

   * - Parameter
     - Type
     - Default
     - Description
   * - ``npix``
     - int
     - ``256``
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
     - list of str
     - ``["xy", "xz", "yz"]``
     - Viewing angles to render. Built-in: ``xy``, ``xz``, ``yz``.
   * - ``custom_views``
     - dict
     - ``{}``
     - User-defined views merged over the built-ins: ``{name: {plane_id, axis1, axis2, theta, phi}}``.
   * - ``fov_au``
     - float or None
     - ``None``
     - [AU] Field of view (full width). None = full grid extent.
   * - ``polaris_cmd``
     - str
     - ``"CMD_DUST_EMISSION"``
     - POLARIS command: ``CMD_DUST_EMISSION``, ``CMD_DUST_SCATTERING``, or both.
   * - ``alignment``
     - str
     - ``"ALIG_PA"``
     - Grain alignment: ``ALIG_PA``, ``ALIG_IDG``, ``ALIG_RAT``, ``ALIG_INTERNAL``.
   * - ``peel_off``
     - bool
     - ``True``
     - Use the peel-off technique for scattering.
   * - ``acceptance_angle``
     - float or None
     - ``None``
     - [deg] Acceptance angle for scattered light. None = POLARIS default.
   * - ``nr_photons_scat``
     - int or None
     - ``None``
     - Photon packages for the scattering Monte-Carlo. None = no scattering source.
   * - ``scat_source_radius_rsun``
     - float
     - ``1.0``
     - [Rsun] Radius of the default scattering source star (used when ``nr_photons_scat`` is set but no explicit source is given).
   * - ``scat_source_temp_k``
     - float
     - ``5000.0``
     - [K] Temperature of the default scattering source star.



Rendering (method arguments)
----------------------------

Passed to ``render_images``. Config-backed arguments default to ``None`` and
resolve from ``cfg.imaging`` (observation), ``cfg.dust`` (dust) and
``cfg.polaris`` (execution). ``custom_views`` and the scattering-source star
parameters are taken directly from ``cfg.imaging``.

.. list-table::
   :header-rows: 1
   :widths: 24 22 16 38

   * - Argument
     - Type
     - Default
     - Description
   * - ``dust_mixtures``
     - dict or None
     - ``None``
     - Dust definitions. If None, uses ``dust.mixtures``.
   * - ``npix``
     - int or None
     - ``None``
     - Image resolution. Falls back to ``imaging.npix``.
   * - ``distance_pc``
     - float or None
     - ``None``
     - Source distance in parsecs. Falls back to ``imaging.distance_pc``.
   * - ``wavelengths_mm``
     - list or None
     - ``None``
     - Wavelengths in mm. Falls back to ``imaging.wavelengths_mm``.
   * - ``views``
     - list or None
     - ``None``
     - Views to render. Falls back to ``imaging.views``.
   * - ``fov_m``
     - float or None
     - ``None``
     - [m] Field of view (full width). Overrides ``fov_au`` if set.
   * - ``fov_au``
     - float or None
     - ``None``
     - [AU] Field of view. Falls back to ``imaging.fov_au``. For sub-boxes, match the RADMC-3D regrid box.
   * - ``label``
     - str
     - ``"whole"``
     - Output subdirectory label (per-call, not config-backed).
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
     - Gas mu. Falls back to ``dust.mean_molecular_weight``.
   * - ``mass_fraction``
     - float or None
     - ``None``
     - Dust-to-gas mass fraction. Falls back to ``dust.mass_fraction``.
   * - ``polaris_binary``
     - str or None
     - ``None``
     - POLARIS executable. Falls back to ``polaris.polaris_binary``.
   * - ``cleanup_views``
     - bool
     - ``True``
     - Remove previous image output for each view before rendering.
   * - ``polaris_cmd``
     - str or None
     - ``None``
     - POLARIS command. Falls back to ``imaging.polaris_cmd``.
   * - ``alignment``
     - str or None
     - ``None``
     - Grain alignment. Falls back to ``imaging.alignment``.
   * - ``peel_off``
     - bool or None
     - ``None``
     - Peel-off technique. Falls back to ``imaging.peel_off``.
   * - ``acceptance_angle``
     - float or None
     - ``None``
     - Acceptance angle for scattered light. Falls back to ``imaging.acceptance_angle``.
   * - ``nr_photons_scat``
     - int or None
     - ``None``
     - Scattering photon packages. Falls back to ``imaging.nr_photons_scat``.
   * - ``source_star_scat``
     - list of dict or None
     - ``None``
     - Explicit stellar sources for scattering. If None and ``nr_photons_scat`` is set, a default source is used (see ``imaging.scat_source_*``).
   * - ``subbox``
     - None / str / list / "all"
     - ``None``
     - Which target to render: None (full cloud), a single sink, a list, or "all" sub-boxes.