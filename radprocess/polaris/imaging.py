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
    n_dust,
    dust_mixtures,
    mean_molecular_weight,
    mass_fraction,
    npix,
    distance_pc,
    wavelengths_mm,
    nr_threads = 1,
    fov_m = None,
    polaris_cmd = "CMD_DUST_EMISSION",
    alignment = "ALIG_PA",
    peel_off = True,
    acceptance_angle = None,
    nr_photons_scat = None,
    source_star_scat = None,
    scat_source_radius_rsun = 1.0,
    scat_source_temp_k = 5000.0,
    detector_shift_m = None,
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
    polaris_cmd : str
        POLARIS command. Options:
        - "CMD_DUST_EMISSION": thermal dust emission (Stokes I, Q, U)
        - "CMD_DUST_SCATTERING": scattered light via Monte Carlo
        - "CMD_DUST_EMISSION CMD_DUST_SCATTERING": both in sequence
    alignment : str
        Grain alignment mechanism. Options:
        - "ALIG_PA": paramagnetic alignment (default)
        - "ALIG_IDG": imperfect Davis-Greenstein
        - "ALIG_RAT": radiative alignment torques
        - "ALIG_INTERNAL": internal alignment
        - "" or None: no alignment (randomly oriented grains)
    peel_off : bool
        If True, use the peel-off technique for scattering. More efficient
        for producing images. Only relevant for CMD_DUST_SCATTERING.
    acceptance_angle : float or None
        Acceptance angle for scattered light in degrees. If None, POLARIS
        uses its default. Only relevant for CMD_DUST_SCATTERING.
    nr_photons_scat : int or None
        Number of photon packages for the scattering Monte Carlo. If None,
        no source_star line is written (only relevant for CMD_DUST_SCATTERING).
    source_star_scat : list of dict or None
        Stellar sources for scattering. Each dict has pos_m, radius_rsun,
        temperature_k. If None and nr_photons_scat is set, a default point
        source at the origin is used.

    Returns
    -------
    cmd_path : Path
    """
    cmd_path = Path(cmd_path)
    cmd_path.parent.mkdir(parents=True, exist_ok=True)

    # size_edges = np.logspace(np.log10(dust_size_min),
    #                          np.log10(dust_size_max), n_dust + 1)

    distance_m = distance_pc * pc2m

    print(f"Writing imaging cmd for '{view_name}' view: {cmd_path}")

    with open(cmd_path, "w") as f:
        # --- <common> block ---
        f.write("<common>\n")

        # for i in range(n_dust):
        #     for comp in dust_components:
        #         f.write(
        #             f'\n\t<dust_component id = "{i}"> '
        #             f'"{comp["path"]}" "plaw" {comp["weight"]} 0 '
        #             f'{size_edges[i]:.2e} {size_edges[i+1]:.2e} '
        #             f'{dust_size_powerlaw}'
        #         )

        for mixture, components in dust_mixtures.items():
            for component in components.values():
                f.write(
                    f'\n\t<dust_component id = "{mixture}"> '
                    f'"{component["path"]}" '
                    f'"{component["distribution"]}" '
                    f'{component["fraction"]} '
                    f'{component["density"]} '
                    f'{component["amin"]:.2e} '
                    f'{component["amax"]:.2e}'
                )
                index = component['index'] if isinstance(component['index'], list) else [component['index']]
                for i in index:
                    f.write(
                        f' {i}'
                    )

        f.write(f"\n\n\t<nr_threads> {nr_threads}\n")
        f.write("\n</common>\n")

        # --- <task> block ---
        f.write("\n<task> 1\n")

        # POLARIS command(s)
        for cmd in polaris_cmd.split():
            f.write(f"\n\t<cmd> {cmd}\n")

        f.write(f'\n\t<path_grid> "{grid_path}"')
        out_str = str(output_path)
        if not out_str.endswith(os.sep):
            out_str += os.sep
        f.write(f'\n\t<path_out>  "{out_str}"\n')

        f.write(f"\n\t<mu> {mean_molecular_weight}\n")
        f.write(f"\n\t<mass_fraction> {mass_fraction}\n")

        # Grain alignment
        if alignment:
            f.write(f"\n\t<align> {alignment}\n")

        # Scattering options
        if peel_off:
            f.write("\n\t<peel_off> 1\n")

        if acceptance_angle is not None:
            f.write(f"\n\t<acc_sca> {acceptance_angle}\n")

        # Scattering source stars
        if nr_photons_scat is not None:
            if source_star_scat is not None:
                for star in source_star_scat:
                    xpos, ypos, zpos = star["pos_m"]
                    r_rsun = star["radius_rsun"]
                    temp = star["temperature_k"]
                    f.write(
                        f'\n\t<source_star nr_photons = "{nr_photons_scat}"> '
                        f"{xpos:17.10e} {ypos:17.10e} {zpos:17.10e} "
                        f"{r_rsun:17.10e} {temp:17.10e}"
                    )
            else:
                # Default point source at origin, with configurable radius/temp.
                f.write(
                    f'\n\t<source_star nr_photons = "{nr_photons_scat}"> '
                    f"0.0 0.0 0.0 {scat_source_radius_rsun} {scat_source_temp_k}"
                )

        plane_id = view_details["plane_id"]
        f.write(f"\n\t<write_inp_midplanes> {npix}")
        f.write(f"\n\t<write_3d_midplanes> {plane_id} {npix}\n")

        # f.write(f"\n\t<midplane_zoom> {midplane_zoom}\n")

        axis1 = view_details["axis1"]
        axis2 = view_details["axis2"]
        f.write(f"\n\t<axis1> {axis1[0]} {axis1[1]} {axis1[2]}")
        f.write(f"\n\t<axis2> {axis2[0]} {axis2[1]} {axis2[2]}\n")

        theta = view_details["theta"]
        phi = view_details["phi"]

        # Compute 2D detector shift from 3D sink offset.
        # POLARIS detector_dust accepts ∆x ∆y as the last two parameters,
        # which shift the image center in the detector plane (in metres).
        det_dx = 0.0
        det_dy = 0.0
        if detector_shift_m is not None:
            shift = np.array(detector_shift_m, dtype=float)
            a1 = np.array(axis1, dtype=float)
            a2 = np.array(axis2, dtype=float)
            det_dx = float(np.dot(shift, a1))
            det_dy = float(np.dot(shift, a2))

        has_shift = (det_dx != 0.0 or det_dy != 0.0)

        for wave_mm in wavelengths_mm:
            wave_m = wave_mm * 1e-3  # mm -> m
            if fov_m is not None:
                line = (
                    f'\n\t<detector_dust nr_pixel = "{npix}"> '
                    f"{wave_m:e} {wave_m:e} 1 1 "
                    f"{theta} {phi} {distance_m:e} {fov_m:e} {fov_m:e}"
                )
                if has_shift:
                    line += f" {det_dx:e} {det_dy:e}"
                f.write(line)
            else:
                if has_shift:
                    # Must provide sidelength to reach the shift fields
                    f.write(
                        f'\n\t<detector_dust nr_pixel = "{npix}"> '
                        f"{wave_m:e} {wave_m:e} 1 1 "
                        f"{theta} {phi} {distance_m:e} -1 -1 {det_dx:e} {det_dy:e}"
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
    n_dust,
    dust_mixtures,
    mean_molecular_weight,
    mass_fraction,
    npix,
    distance_pc,
    wavelengths_mm,
    views = None,
    custom_views = None,
    nr_threads = 1,
    fov_m = None,
    output_num = 0,
    polaris_binary = "polaris",
    label = "whole",
    cleanup_views = True,
    polaris_cmd = "CMD_DUST_EMISSION",
    alignment = "ALIG_PA",
    peel_off = True,
    acceptance_angle = None,
    nr_photons_scat = None,
    source_star_scat = None,
    scat_source_radius_rsun = 1.0,
    scat_source_temp_k = 5000.0,
    detector_shift_m = None,
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

    # User-defined views override / extend the built-in axis-aligned ones.
    available_views = dict(STANDARD_VIEWS)
    if custom_views:
        available_views.update(custom_views)

    print(f"\nRendering {len(views)} views: {views}")
    print(f"  Wavelengths: {wavelengths_mm} mm")
    print(f"  Resolution:  {npix} x {npix} px")
    print(f"  Distance:    {distance_pc} pc")
    print(f"  Label:       {label}\n")

    image_dirs = {}

    for view_name in views:
        if view_name not in available_views:
            raise ValueError(
                f"Unknown view '{view_name}'. "
                f"Available: {list(available_views.keys())}"
            )

        view_details = available_views[view_name]
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
            cmd_path = cmd_path,
            grid_path = grid_path,
            output_path = view_output,
            view_name = view_name,
            view_details = view_details,
            n_dust = n_dust,
            dust_mixtures = dust_mixtures,
            mean_molecular_weight = mean_molecular_weight,
            mass_fraction = mass_fraction,
            npix = npix,
            distance_pc = distance_pc,
            wavelengths_mm = wavelengths_mm,
            nr_threads = nr_threads,
            fov_m = fov_m,
            polaris_cmd = polaris_cmd,
            alignment = alignment,
            peel_off = peel_off,
            acceptance_angle = acceptance_angle,
            nr_photons_scat = nr_photons_scat,
            source_star_scat = source_star_scat,
            scat_source_radius_rsun = scat_source_radius_rsun,
            scat_source_temp_k = scat_source_temp_k,
            detector_shift_m = detector_shift_m,
        )

        # Run POLARIS
        log_filename = f"polaris_render_{view_name}_{output_num:05d}.{label}.log"
        log_path = polaris_dir / log_filename

        print(f"\nRendering '{view_name}' ({label})...")
        run_polaris(cmd_path, log_path=log_path, polaris_binary=polaris_binary)

        image_dirs[view_name] = view_output

    print(f"All {label} images rendered successfully!")
    print(f"Output directories: {image_dirs}\n")

    return image_dirs