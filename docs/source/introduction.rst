.. _chap-introduction:

Introduction
************

What is radprocess?
===================

Radprocess is a Python package for post-processing magnetohydrodynamic (MHD)
simulations performed with the `RAMSES <https://bitbucket.org/rteyssie/ramses>`_
adaptive mesh refinement code. It provides a unified pipeline to convert RAMSES
outputs into input files for two radiative transfer codes,
`RADMC-3D <https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/>`_
and `POLARIS <https://portia.astrophysik.uni-kiel.de/polaris/>`_, and to
orchestrate the full chain from raw simulation data to synthetic dust continuum
images.

The package was developed in the context of the
`ERC Synergy ECOGAL <https://www.ecogal.eu/>`_ project, with a focus on
protostellar collapse and protoplanetary disk formation in multi-fluid
(gas + dust) RAMSES simulations.


Why radprocess?
===============

Producing synthetic observations from MHD simulations typically requires a long
sequence of format conversions, unit transformations, opacity computations, and
code executions that are error-prone when done by hand. Radprocess automates
this chain into a reproducible, step-by-step pipeline.

The pipeline supports simulations with an arbitrary number of dust species
(size bins), handles the octree AMR structure natively (no regridding to a
uniform mesh), and takes care of the CGS/SI unit conversions required by
RADMC-3D and POLARIS respectively.


The pipeline at a glance
========================

The full pipeline consists of eight steps:

1. **Load RAMSES** -- Read the RAMSES AMR output via pymses and store the
   cell data (density, velocity, temperature, magnetic field) in a compressed
   Zarr archive for efficient reuse.

2. **Convert to POLARIS grid** -- Build a POLARIS-format binary octree with densities 
   converted to SI units (kg/m\ :sup:`3`).

3. **Convert to RADMC-3D grid** -- Build a RADMC-3D-format octree
   (``amr_grid.inp``, ``dust_density.inp``, ``stars.inp``), keeping CGS units (g/cm\ :sup:`3`).

4. **Run POLARIS opacity** -- Execute POLARIS with a single photon package
   to generate wavelength-dependent dust opacity tables
   (``dust_mixture_*.dat``) for each size bin.

5. **Prepare RADMC-3D inputs** -- Convert the POLARIS opacity tables to
   RADMC-3D format (``dustkappa_*.inp``) and write the remaining input files
   (``dustopac.inp``, ``wavelength_micron.inp``, ``radmc3d.inp``).

6. **Run RADMC-3D mctherm** -- Execute the RADMC-3D Monte Carlo thermal
   computation to obtain the dust temperature field
   (``dust_temperature.bdat``).

7. **Merge temperature** -- Inject the RADMC-3D dust temperatures back into
   the POLARIS grid, producing a merged grid (``grid_temp.radmc3d.dat``) that
   is physically consistent with the RADMC-3D thermal solution.

8. **Render images** -- Run POLARIS ``CMD_DUST_EMISSION`` on the merged grid
   to produce synthetic dust continuum images at user-specified wavelengths
   and viewing angles.

Each step can be run independently. The user does not need to start from scratch at each session.
There also exist additional actions (e.g. rending images with RADMC3D, post-processing inside a subbox of the AMR grid, etc.)
that are described later in the documentation.



Acknowledgements
================

If you use radprocess in a publication, please cite the relevant papers
(citation information will be provided here once available).


Copyright
=========

Radprocess is distributed under the MIT License.
It was originally developed for the ENYGMA Large Program (NOEMA) under the ERC Synergy Grant
“ECOGAL” (project ID 855130). 