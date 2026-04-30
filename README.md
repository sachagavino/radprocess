# radprocess Version 1.0.0

## Introduction
RadProcess is a pipeline that converts RAMSES output for radiative tranfer postprocessing. 
It also uses the web interface Gradio, which provides an intuitive way to build and run your postprocessing pipeline, but it is usable from a notebook.
RadProcess uses a combination between Zarr and Xarray for a fast and easy way to handle the data.

Here is the pipeline dashboard:

<picture>
  <source srcset="radprocess/interface/_static/dashboad_gr_pipeline.png" media="(prefers-color-scheme: dark)">
  <img src="radprocess/interface/_static/dashboad_gr_pipeline.png" alt="dashboard", width="700">
</picture>

## User's guide
The file INSTALL.md provides the full installation procedure. Alternatively, it is also possible to read the installation procedure and a how-to-use Radprocess in the documentation [here][1].

## External requirements (not installable via pip):

- RAMSES
- RADMC3D
- POLARIS
- PyMSES (see custom fork)

## Branches
main: stable up-to-date version.

dev: use it at your own risks.

## News
[11.06.2025]: Includes multiple dust populations


[1]: https://radprocess.readthedocs.io/en/latest
