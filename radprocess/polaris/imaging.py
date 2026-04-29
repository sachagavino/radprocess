"""
_____________________________________________________________________________________________________________
file name: imaging
last update: Apr 2026
language: > PYTHON 3.9
short description: Generate POLARIS command files for dust continuum imaging (CMD_DUST_EMISSION)
                   and execute POLARIS to produce synthetic observations at specified
                   wavelengths and viewing angles.

                   This is Step 8 of the pipeline. It uses the merged POLARIS grid
                   (grid_temp.radmc3d.dat from Step 7) which carries the RADMC-3D thermal
                   solution, ensuring physical consistency between the temperature
                   computation and the final images.
_____________________________________________________________________________________________________________
"""

import os
import sys
import shutil
import subprocess
from pathlib import Path

import numpy as np

from radprocess.constants.constants import pc2m
from radprocess.polaris.opacity import _find_info_file, _read_ndust, run_polaris


# ============================================================
#  Standard viewing geometries
# ============================================================

# Each view is defined by:
#   plane_id : int used by POLARIS for write_3d_midplanes
#   axis1, axis2 : detector orientation vectors
#   theta, phi : observer angles in degrees

STANDARD_VIEWS = {
    "xy": {
        "plane_id": 1,
        "axis1": [1, 0, 0],
        "axis2": [0, 1, 0],
        "theta": 0,
        "phi": 0,
    },
    "xz": {
        "plane_id": 2,
        "axis1": [1, 0, 0],
        "axis2": [0, 0, 1],
        "theta": 90,
        "phi": 0,
    },
    "yz": {
        "plane_id": 3,
        "axis1": [0, 1, 0],
        "axis2": [0, 0, 1],
        "theta": 90,
        "phi": 90,
    },
}


# ============================================================
#  Write imaging command file
# ============================================================

def write_imaging_cmd(
    cmd_path,
    grid_path,
    output_path,
    view_name,
    view_details,
    dust_components,
    n_dust,
    dust_size_min,
    dust_size_max,
    dust_size_powerlaw,
    mean_molecular_weight,
    mass_fraction,
    npix,
    distance_pc,
    wavelengths_mm,
    nr_threads=8,
    midplane_zoom=1,
    fov_m=None,
):
    """
    Write a POLARIS command file for dust continuum imaging at a
    specific viewing angle.

    Parameters
    ----------
    cmd_path : str or Path
        Where to write the .cmd file.
    grid_path : str or Path
        Path to the merged POLARIS grid (grid_temp.radmc3d.dat).
    output_path : str or Path
        Directory where POLARIS will write image output.
    view_name : str
        Name of the view (e.g. "xy", "xz", "yz").
    view_details : dict
        Viewing geometry with keys: plane_id, axis1, axis2, theta, phi.
    dust_components : list of dict
        Dust material definitions (same as for opacity run).
    n_dust : int
        Number of dust species.
    dust_size_min, dust_size_max : float
        Grain size range in metres.
    dust_size_powerlaw : float
        Size distribution power-law exponent.
    mean_molecular_weight : float
        Gas mu.
    mass_fraction : float
        Dust-to-gas mass fraction.
    npix : int
        Image resolution (npix x npix).
    distance_pc : float
        Source distance in parsecs.
    wavelengths_mm : list of float
        Wavelengths to image, in millimetres.
    nr_threads : int
        Number of OpenMP threads.
    midplane_zoom : int or float
        Zoom factor for midplane output (1 = full box, 10 = zoomed).
    fov_m : float or None
        Field of view in metres (half-width). If None, POLARIS uses
        the full grid extent (for whole-box imaging). If provided,
        the detector_dust line includes explicit sidelength parameters
        (for zoomed/inner imaging).

    Returns
    -------
    cmd_path : Path
    """
    cmd_path = Path(cmd_path)
    cmd_path.parent.mkdir(parents=True, exist_ok=True)

    size_edges = np.logspace(np.log10(dust_size_min),
                             np.log10(dust_size_max), n_dust + 1)

    distance_m = distance_pc * pc2m

    print(f"  Writing imaging cmd for '{view_name}' view: {cmd_path}")

    with open(cmd_path, "w") as f:
        # --- <common> block ---
        f.write("<common>\n")

        for i in range(n_dust):
            for comp in dust_components:
                f.write(
                    f'\n\t<dust_component id = "{i}"> '
                    f'"{comp["path"]}" "plaw" {comp["weight"]} 0 '
                    f'{size_edges[i]:.2e} {size_edges[i+1]:.2e} '
                    f'{dust_size_powerlaw}'
                )

        f.write(f"\n\n\t<nr_threads> {nr_threads}\n")
        f.write("\n</common>\n")

        # --- <task> block ---
        f.write("\n<task> 1\n")
        f.write("\n\t<cmd> CMD_DUST_EMISSION\n")

        f.write(f'\n\t<path_grid> "{grid_path}"')
        out_str = str(output_path)
        if not out_str.endswith(os.sep):
            out_str += os.sep
        f.write(f'\n\t<path_out>  "{out_str}"\n')

        f.write(f"\n\t<mu> {mean_molecular_weight}\n")
        f.write(f"\n\t<mass_fraction> {mass_fraction}\n")

        f.write("\n\t<align> ALIG_PA\n")

        plane_id = view_details["plane_id"]
        f.write(f"\n\t<write_inp_midplanes> {npix}")
        f.write(f"\n\t<write_3d_midplanes> {plane_id} {npix}\n")

        f.write(f"\n\t<midplane_zoom> {midplane_zoom}\n")

        axis1 = view_details["axis1"]
        axis2 = view_details["axis2"]
        f.write(f"\n\t<axis1> {axis1[0]} {axis1[1]} {axis1[2]}")
        f.write(f"\n\t<axis2> {axis2[0]} {axis2[1]} {axis2[2]}\n")

        theta = view_details["theta"]
        phi = view_details["phi"]

        for wave_mm in wavelengths_mm:
            wave_m = wave_mm * 1e-3  # mm -> m
            if fov_m is not None:
                f.write(
                    f'\n\t<detector_dust nr_pixel = "{npix}"> '
                    f"{wave_m:e} {wave_m:e} 1 1 "
                    f"{theta} {phi} {distance_m:e} {fov_m:e} {fov_m:e}"
                )
            else:
                f.write(
                    f'\n\t<detector_dust nr_pixel = "{npix}"> '
                    f"{wave_m:e} {wave_m:e} 1 1 "
                    f"{theta} {phi} {distance_m:e}"
                )

        f.write("\n\n</task>")

    return cmd_path


# ============================================================
#  Render images for one or more views
# ============================================================

def render_images(
    polaris_dir,
    image_output_dir,
    grid_path,
    dust_components,
    n_dust,
    dust_size_min,
    dust_size_max,
    dust_size_powerlaw,
    mean_molecular_weight,
    mass_fraction,
    npix,
    distance_pc,
    wavelengths_mm,
    views=None,
    nr_threads=8,
    midplane_zoom=1,
    fov_m=None,
    output_num=0,
    polaris_binary="polaris",
    label="whole",
    cleanup_views=True,
):
    """
    Full Step 8: write POLARIS imaging command files and execute them
    for each requested viewing angle.

    Parameters
    ----------
    polaris_dir : str or Path
        POLARIS working directory (for storing .cmd and .log files).
    image_output_dir : str or Path
        Base directory for image output. Subdirectories are created per
        view (e.g. {image_output_dir}/{label}/xy/).
    grid_path : str or Path
        Path to the merged POLARIS grid (grid_temp.radmc3d.dat).
    dust_components : list of dict
        Dust material definitions.
    n_dust : int
        Number of dust species.
    dust_size_min, dust_size_max : float
        Grain size range in metres.
    dust_size_powerlaw : float
        Size distribution power-law exponent.
    mean_molecular_weight : float
        Gas mu.
    mass_fraction : float
        Dust-to-gas mass fraction.
    npix : int
        Image resolution (npix x npix).
    distance_pc : float
        Source distance in parsecs.
    wavelengths_mm : list of float
        Wavelengths to image in mm.
    views : list of str or None
        Which views to render (e.g. ["xy", "xz", "yz"]).
        If None, renders all three standard views.
    nr_threads : int
        Number of OpenMP threads.
    midplane_zoom : int or float
        Midplane zoom factor.
    fov_m : float or None
        Field of view in metres. None = full box.
    output_num : int
        RAMSES output number (used for naming).
    polaris_binary : str
        Name or path of the POLARIS executable.
    label : str
        Subdirectory label ("whole" or "inner").
    cleanup_views : bool
        If True, remove previous image output for each view before rendering.

    Returns
    -------
    image_dirs : dict
        {view_name: Path} mapping each view to its output directory.
    """
    polaris_dir = Path(polaris_dir)
    image_output_dir = Path(image_output_dir)

    if views is None:
        views = ["xy", "xz", "yz"]

    print(f"\nRendering {len(views)} views: {views}")
    print(f"  Wavelengths: {wavelengths_mm} mm")
    print(f"  Resolution:  {npix} x {npix} px")
    print(f"  Distance:    {distance_pc} pc")
    print(f"  Zoom:        {midplane_zoom}x")
    print(f"  Label:       {label}\n")

    image_dirs = {}

    for view_name in views:
        if view_name not in STANDARD_VIEWS:
            raise ValueError(
                f"Unknown view '{view_name}'. "
                f"Available: {list(STANDARD_VIEWS.keys())}"
            )

        view_details = STANDARD_VIEWS[view_name]
        view_output = image_output_dir / label / view_name
        view_output.mkdir(parents=True, exist_ok=True)

        # Cleanup previous output
        if cleanup_views and view_output.exists():
            shutil.rmtree(view_output)
            view_output.mkdir(parents=True, exist_ok=True)

        # Write cmd file
        cmd_filename = f"polaris_render_{view_name}_{output_num:05d}.{label}.cmd"
        cmd_path = polaris_dir / cmd_filename

        write_imaging_cmd(
            cmd_path=cmd_path,
            grid_path=grid_path,
            output_path=view_output,
            view_name=view_name,
            view_details=view_details,
            dust_components=dust_components,
            n_dust=n_dust,
            dust_size_min=dust_size_min,
            dust_size_max=dust_size_max,
            dust_size_powerlaw=dust_size_powerlaw,
            mean_molecular_weight=mean_molecular_weight,
            mass_fraction=mass_fraction,
            npix=npix,
            distance_pc=distance_pc,
            wavelengths_mm=wavelengths_mm,
            nr_threads=nr_threads,
            midplane_zoom=midplane_zoom,
            fov_m=fov_m,
        )

        # Run POLARIS
        log_filename = f"polaris_render_{view_name}_{output_num:05d}.{label}.log"
        log_path = polaris_dir / log_filename

        print(f"\n  Rendering '{view_name}' ({label})...")
        run_polaris(cmd_path, log_path=log_path, polaris_binary=polaris_binary)

        image_dirs[view_name] = view_output

    print(f"\nAll {label} images rendered successfully!")
    print(f"Output directories: {image_dirs}\n")

    return image_dirs