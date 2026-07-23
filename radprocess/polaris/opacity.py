"""
_____________________________________________________________________________________________________________
file name: opacity
last update: Apr 2026
language: > PYTHON 3.9
short description: Generate a POLARIS command file for the opacity-only run (CMD_TEMP with
                   1 photon package per star) and execute POLARIS. The sole output of this
                   step is the dust opacity table (dust_mixture_*.dat) that will later be
                   converted into RADMC-3D format.
_____________________________________________________________________________________________________________
"""

import os
import sys
import shutil
import subprocess
from pathlib import Path

import numpy as np

from radprocess.constants.constants import (
    Ggram, M_sun, L_sun, R_sun, sigma
)
from radprocess import ramses


# ============================================================
#  Derive stellar properties from RAMSES sinks
# ============================================================

def derive_stars_properties(ramses_dir, f_acc=0.1):
    """
    Read the RAMSES sink file and derive stellar properties
    (position, mass, radius, luminosity, Teff) for every sink
    that has a nonzero internal luminosity.

    Parameters
    ----------
    ramses_dir : str or Path
        Path to the RAMSES output directory.
    f_acc : float
        Accretion efficiency factor.

    Returns
    -------
    stars : list of dict
        Each dict contains: id, mass_msun, pos_m (list of 3 floats in metres,
        centred on the domain), radius_rsun, luminosity_lsun, temperature_k.
    boxlen : float
        Domain box length in parsecs (used by callers for coordinate transforms).
    """
    sinks = ramses.read.sink_info(ramses_dir)

    if sinks.num_sinks == 0:
        return [], 0.0

    cols = sinks.columns
    id_col = cols.index("Id")
    x_col = cols.index("x")
    y_col = cols.index("y")
    z_col = cols.index("z")
    m_col = cols.index("M[Msol]")
    accrate_col = cols.index("acc_rate[Msol/y]")
    acclum_col = cols.index("acc_lum[Lsol]")
    intlum_col = cols.index("int_lum[Lsol]")
    teff_col = cols.index("Teff[K]")

    sec_yr = 365.25 * 24 * 3600
    data = sinks.data

    # Infer boxlen from the Zarr or from the domain size
    # For now, compute from sink positions (they are in code units = pc)
    # We need to read it from the info file instead
    info_path = _find_info_file(ramses_dir)
    boxlen = _read_boxlen(info_path)
    unit_l_cm = _read_unit_l(info_path)
    unit_l_m  = unit_l_cm * 1e-2

    stars = []
    for i in range(sinks.num_sinks):
        lint = data[i, intlum_col]
        lacc = data[i, acclum_col]
        teff = data[i, teff_col]

        # Skip sinks without luminosity info
        if lint == 0 and lacc == 0:
            continue

        ltot = (lint + lacc) * L_sun  # erg/s
        teff_K = teff  # K

        mass_g = data[i, m_col] * M_sun
        acc_rate_g = data[i, accrate_col] * M_sun / sec_yr
        lacc_cgs = lacc * L_sun

        if lacc_cgs > 0 and teff_K > 0:
            radius_cm = np.sqrt(ltot / (4 * np.pi * sigma * teff_K**4))
        else:
            radius_cm = 1.0 * R_sun


        # Position centred on the domain, in metres
        pos_code = np.array([data[i, x_col], data[i, y_col], data[i, z_col]])
        pos_m = (pos_code - boxlen / 2.0) * unit_l_m

        #pos_m = (pos_pc - boxlen_pc / 2.0) * pc2m

        stars.append({
            "id": int(data[i, id_col]),
            "mass_msun": data[i, m_col],
            "pos_m": pos_m.tolist(),
            "radius_rsun": radius_cm / R_sun,
            "luminosity_lsun": (lint + lacc),
            "temperature_k": teff_K,
        })

    print(f"Derived properties for {len(stars)} stars "
          f"(skipped {sinks.num_sinks - len(stars)} without luminosity).")

    return stars, boxlen


# ============================================================
#  Internal helpers
# ============================================================

def _find_info_file(ramses_dir):
    """Locate info_XXXXX.txt in a RAMSES output directory."""
    ramses_dir = Path(ramses_dir)
    # Extract output number from directory name
    dirname = ramses_dir.name  # e.g. "output_00940"
    num_str = dirname.split("_")[-1]
    info_file = ramses_dir / f"info_{num_str}.txt"
    if info_file.exists():
        return info_file
    raise FileNotFoundError(f"Info file not found: {info_file}")


def _read_boxlen(info_path):
    """Read boxlen from a RAMSES info file."""
    with open(info_path, "r") as f:
        for line in f:
            if "boxlen" in line and "=" in line:
                return float(line.split("=")[1])
    raise ValueError(f"'boxlen' not found in {info_path}")


def _read_ndust(info_path):
    """Read ndust from a RAMSES info file. Returns 0 if not found."""
    with open(info_path, "r") as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("ndust"):
                return int(stripped.split("=")[1])
    return 0

def _read_unit_l(info_path):
    """Read unit_l (length conversion factor, code units → cm) from a RAMSES info file."""
    with open(info_path, "r") as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("unit_l") and "=" in stripped:
                return float(stripped.split("=")[1])
    raise ValueError(f"'unit_l' not found in {info_path}")


# ============================================================
#  Cleanup
# ============================================================

def cleanup_previous_run(output_path):
    """Remove output files and directories from a previous POLARIS run."""
    output_path = Path(output_path)

    for subdir in ["data", "plots"]:
        d = output_path / subdir
        if d.exists():
            shutil.rmtree(d)
            print(f"  Removed directory: {d}")

    grid_temp = output_path / "grid_temp.dat"
    if grid_temp.exists():
        grid_temp.unlink()
        print(f"  Removed file: {grid_temp}")


# ============================================================
#  Write POLARIS .cmd file
# ============================================================

def write_opacity_cmd(
    cmd_path,
    grid_path,
    output_path,
    stars,
    n_dust,
    # dust_components,
    # dust_size_min,
    # dust_size_max,
    # dust_size_powerlaw = -3.5,
    dust_mixtures,
    mean_molecular_weight = 2.37,
    mass_fraction = 1,
    nr_threads = 1,
    output_num = 0,
):
    """
    Write a POLARIS command file for the opacity-only run.

    Parameters
    ----------
    cmd_path : str or Path
        Where to write the .cmd file.
    grid_path : str or Path
        Path to the POLARIS binary octree grid file.
    output_path : str or Path
        Directory where POLARIS will write its output (data/ subdirectory).
    stars : list of dict
        Stellar properties as returned by derive_stars_properties().
    n_dust : int
        Number of dust species (size bins).
    dust_components : list of dict
        Each dict describes a dust material component with keys:
            path (str): path to the .nk or .cs file
            weight (float): mass fraction weight (e.g. 0.625 for silicate)
        Example for silicate + carbon:
            [{"path": "/path/silicate.cs", "weight": 0.625},
             {"path": "/path/carbon.cs",   "weight": 0.375}]
        For a single component, just one entry with weight=1.0.
    dust_size_min : float
        Minimum grain radius in metres.
    dust_size_max : float
        Maximum grain radius in metres.
    dust_size_powerlaw : float
        Power-law exponent for the size distribution (e.g. -3.5 for MRN).
    mean_molecular_weight : float
        Mean molecular weight (mu).
    mass_fraction : float
        Dust-to-gas mass fraction.
    nr_threads : int
        Number of OpenMP threads for POLARIS.
    output_num : int
        RAMSES output number (used only for naming).
    """
    cmd_path = Path(cmd_path)
    cmd_path.parent.mkdir(parents=True, exist_ok=True)

    # Log-spaced size bin edges
    # size_edges = np.logspace(np.log10(dust_size_min),
    #                          np.log10(dust_size_max), n_dust + 1)

    print(f"Writing POLARIS opacity command file: {cmd_path}")

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
        f.write("\n\t<cmd> CMD_TEMP\n")

        for star in stars:
            xpos, ypos, zpos = star["pos_m"]
            r_rsun = star["radius_rsun"]
            temp = star["temperature_k"]
            f.write(
                f"\n\t<source_star nr_photons = \"1\"> "
                f"{xpos:17.10e} {ypos:17.10e} {zpos:17.10e} "
                f"{r_rsun:17.10e} {temp:17.10e}"
            )

        # f.write('\n\t<source_isrf nr_photons = "1">  1\n')

        f.write(f'\n\t<path_grid> "{grid_path}"')
        # POLARIS expects trailing slash on output path
        out_str = str(output_path)
        if not out_str.endswith(os.sep):
            out_str += os.sep
        f.write(f'\n\t<path_out>  "{out_str}"\n')

        f.write(f"\n\t<mu> {mean_molecular_weight}\n")
        f.write(f"\n\t<mass_fraction> {mass_fraction}\n")
        f.write("\n</task>")

    print(f"Command file written: {cmd_path}\n")
    return cmd_path


# ============================================================
#  Run POLARIS
# ============================================================

def run_polaris(cmd_path, log_path=None, polaris_binary="polaris"):
    """
    Execute POLARIS with the given command file.

    Parameters
    ----------
    cmd_path : str or Path
        Path to the .cmd file.
    log_path : str or Path or None
        If provided, POLARIS stdout/stderr is also written to this file.
    polaris_binary : str
        Name or path of the POLARIS executable.

    Raises
    ------
    subprocess.CalledProcessError
        If POLARIS exits with a nonzero return code.
    """
    cmd_path = Path(cmd_path)
    command = [polaris_binary, str(cmd_path)]

    print(f"Executing: {' '.join(command)}")
    if log_path:
        print(f"Log file: {log_path}\n")

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
    )

    log_file = open(log_path, "w") if log_path else None

    try:
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            if log_file:
                log_file.write(line.rstrip("\r\n") + "\n")
    finally:
        if log_file:
            log_file.close()

    process.wait()
    if process.returncode != 0:
        raise subprocess.CalledProcessError(process.returncode, command)


# ============================================================
#  High-level: write cmd + run
# ============================================================

def run_opacity(
    ramses_dir,
    polaris_dir,
    grid_path,
    # dust_components,
    # dust_size_min,
    # dust_size_max,
    # dust_size_powerlaw = -3.5,
    dust_mixtures,
    mean_molecular_weight = 2.37,
    mass_fraction = 0.01,
    nr_threads = 1,
    f_acc = 0.1,
    n_dust_override = None,
    polaris_binary = "polaris",
    cleanup = True,
):
    """
    Full Step 4: derive stellar properties, write the POLARIS command file,
    clean up any previous run, and execute POLARIS.

    Parameters
    ----------
    ramses_dir : str or Path
        RAMSES output directory.
    polaris_dir : str or Path
        POLARIS working directory (where the grid file lives and where
        the .cmd, log, and data/ output will go).
    grid_path : str or Path
        Path to the POLARIS binary octree grid file.
    dust_components : list of dict
        Dust material definitions (see write_opacity_cmd docstring).
    dust_size_min, dust_size_max : float
        Grain size range in metres.
    dust_size_powerlaw : float
        MRN exponent (default -3.5).
    mean_molecular_weight : float
        Gas mu.
    mass_fraction : float
        Dust-to-gas mass fraction.
    nr_threads : int
        OpenMP threads for POLARIS.
    f_acc : float
        Accretion efficiency for deriving stellar radii.
    n_dust_override : int or None
        If provided, use this as the number of dust species instead of
        reading from the RAMSES info file.
    polaris_binary : str
        Name or path to the POLARIS executable.
    cleanup : bool
        If True, remove data/, plots/, grid_temp.dat from a previous run.

    Returns
    -------
    data_dir : Path
        Path to the POLARIS output data/ directory containing
        dust_mixture_*.dat files.
    """
    ramses_dir = Path(ramses_dir)
    polaris_dir = Path(polaris_dir)
    grid_path = Path(grid_path)

    # 1) Derive stellar properties
    stars, boxlen_pc = derive_stars_properties(ramses_dir, f_acc=f_acc)
    if not stars:
        raise ValueError("No stars with luminosity data found in RAMSES sinks.")

    # 2) Determine n_dust
    if n_dust_override is not None:
        n_dust = n_dust_override
    else:
        info_path = _find_info_file(ramses_dir)
        n_dust = _read_ndust(info_path)
        if n_dust == 0:
            n_dust = 1
            print("No dust in simulation. Using a single virtual dust species.")

    # 3) Infer output number
    dirname = ramses_dir.name
    output_num = int(dirname.split("_")[-1])

    # 4) Write .cmd file
    cmd_filename = f"polaris_opacity_{output_num:05d}.cmd"
    cmd_path = polaris_dir / cmd_filename

    write_opacity_cmd(
        cmd_path = cmd_path,
        grid_path = grid_path,
        output_path = polaris_dir,
        stars = stars,
        n_dust = n_dust,
        # dust_components = dust_components,
        # dust_size_min = dust_size_min,
        # dust_size_max = dust_size_max,
        # dust_size_powerlaw = dust_size_powerlaw,
        dust_mixtures = dust_mixtures,
        mean_molecular_weight = mean_molecular_weight,
        mass_fraction = mass_fraction,
        nr_threads = nr_threads,
        output_num = output_num
    )

    # 5) Cleanup previous run
    if cleanup:
        cleanup_previous_run(polaris_dir)

    # 6) Run POLARIS
    log_path = polaris_dir / f"polaris_opacity_{output_num:05d}.log"
    run_polaris(cmd_path, log_path=log_path, polaris_binary=polaris_binary)

    data_dir = polaris_dir / "data"
    return data_dir
