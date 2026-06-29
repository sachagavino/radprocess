"""
_____________________________________________________________________________________________________________
file name: prepare
last update: Apr 2026
language: > PYTHON 3.9
short description: Prepare all remaining RADMC-3D input files after the grid conversion
                   (Step 3) and the POLARIS opacity run (Step 4).

                   This module converts POLARIS dust_mixture_*.dat opacity tables into
                   RADMC-3D dustkappa_*.inp format, and writes dustopac.inp,
                   wavelength_micron.inp, stars.inp, and radmc3d.inp.
_____________________________________________________________________________________________________________
"""

import os
from pathlib import Path

import numpy as np

from radprocess.constants.constants import M_sun, R_sun, pc2cm
from radprocess.polaris.opacity import derive_stars_properties


# ============================================================
#  Convert POLARIS opacities to RADMC-3D format
# ============================================================

def convert_polaris_opacities(polaris_data_dir, radmc_dir, n_dust=None):
    """
    Convert POLARIS dust_mixture_*.dat opacity files into RADMC-3D
    dustkappa_*.inp files.

    Auto-detects the POLARIS output filenames, which may follow either:
        - Old convention: dust_mixture_001.dat, dust_mixture_002.dat, ...
        - New convention: dust_mixture_001_comp_001.dat, ...

    Header lines are auto-detected (lines starting with '#' or that
    cannot be parsed as numbers are skipped).

    The POLARIS opacity table columns (after the header) are:
        wavelength [m], ..., kabs_x, kabs_y, ksca_x, ksca_y

    The last 4 columns are absorption and scattering opacities for
    two polarisation directions. We average them to get orientation-
    averaged values, then convert from POLARIS units (m^2/kg) to
    RADMC-3D units (cm^2/g) by multiplying by 10.

    Parameters
    ----------
    polaris_data_dir : str or Path
        Path to the POLARIS data/ directory containing dust_mixture_*.dat.
    radmc_dir : str or Path
        Output directory where dustkappa_*.inp files will be written.
    n_dust : int or None
        Number of dust species. If None, auto-detected from files found.

    Returns
    -------
    kappa_model_names : list of str
        Model names used in dustopac.inp.
    """
    polaris_data_dir = Path(polaris_data_dir)
    radmc_dir = Path(radmc_dir)
    radmc_dir.mkdir(parents=True, exist_ok=True)

    # Auto-detect POLARIS opacity files
    polaris_files = sorted(polaris_data_dir.glob("dust_mixture_*.dat"))

    if not polaris_files:
        raise FileNotFoundError(
            f"No dust_mixture_*.dat files found in {polaris_data_dir}. "
            "Run the POLARIS opacity step first (Step 4)."
        )

    if n_dust is not None and len(polaris_files) != n_dust:
        print(f"WARNING: Expected {n_dust} dust species but found "
              f"{len(polaris_files)} opacity files. Using found files.")

    kappa_model_names = []

    print("Converting POLARIS opacities to RADMC-3D format:")

    for polaris_file in polaris_files:
        # Derive model name from filename (strip .dat)
        model_name = polaris_file.stem  # e.g. "dust_mixture_001_comp_001"
        kappa_filename = f"dustkappa_{model_name}.inp"
        kappa_model_names.append(model_name)

        # Read data, skipping any line that is a comment (#), empty,
        # or has fewer columns than the actual opacity table.
        rows = []
        n_header = 0
        with open(polaris_file, "r") as pf:
            for line in pf:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    n_header += 1
                    continue
                parts = stripped.split()
                # Opacity data rows have many columns (≥5);
                # skip stray metadata lines (1-2 columns)
                if len(parts) < 5:
                    n_header += 1
                    continue
                try:
                    rows.append([float(p) for p in parts])
                except ValueError:
                    n_header += 1
                    continue

        data = np.array(rows)
        print(f"    {polaris_file.name}: skipped {n_header} header lines, "
              f"{len(data)} wavelength points, {data.shape[1]} columns")

        output_path = radmc_dir / kappa_filename
        with open(output_path, "w") as f:
            # Format 2: lambda, kappa_abs, kappa_sca (no g)
            f.write("2\n")
            f.write(f"{len(data)}\n")

            ncols = data.shape[1]

            for row in data:
                # Wavelength: m -> micron
                wave_um = row[0] / 1e-6

                if ncols >= 21:
                    # New POLARIS format (21 columns):
                    # cols 17,18 = avgKabs1, avgKabs2 [m^2/kg]
                    # cols 19,20 = avgKsca1, avgKsca2 [m^2/kg]
                    # Already mass opacities, convert m^2/kg -> cm^2/g (*10)
                    kabs = (row[17] * 2 + row[18]) / 3 * 10.0
                    ksca = (row[19] * 2 + row[20]) / 3 * 10.0
                else:
                    # Old POLARIS format (fewer columns):
                    # Last 4 cols = kabs_x, kabs_y, ksca_x, ksca_y [m^2/kg]
                    kabs = (row[-4] * 2 + row[-3]) / 3 * 10.0
                    ksca = (row[-2] * 2 + row[-1]) / 3 * 10.0

                f.write(f"{wave_um:e} {kabs:e} {ksca:e}\n")

        print(f"    {kappa_filename}")

    return kappa_model_names


# ============================================================
#  Write dustopac.inp
# ============================================================

def write_dustopac(radmc_dir, kappa_model_names):
    """
    Write the RADMC-3D dustopac.inp index file.

    Parameters
    ----------
    radmc_dir : str or Path
        Directory where dustopac.inp will be written.
    kappa_model_names : list of str
        Model names matching dustkappa_<name>.inp files.
    """
    radmc_dir = Path(radmc_dir)
    filepath = radmc_dir / "dustopac.inp"

    with open(filepath, "w") as f:
        f.write("2\n")
        f.write(f"{len(kappa_model_names)}\n")
        for name in kappa_model_names:
            f.write("-\n")
            f.write("1\n")    # input style: dustkappa_*.inp
            f.write("0\n")    # normal thermal grains
            f.write(f"{name}\n")

    print(f"    dustopac.inp ({len(kappa_model_names)} species)")


# ============================================================
#  Write wavelength_micron.inp
# ============================================================

def write_wavelength_micron(radmc_dir, wave_min=0.27, wave_max=3000, n_wave=200):
    """
    Write the RADMC-3D wavelength grid file.

    Parameters
    ----------
    radmc_dir : str or Path
        Output directory.
    wave_min, wave_max : float
        Wavelength range in microns.
    n_wave : int
        Number of log-spaced wavelength points.

    Returns
    -------
    wavelengths : np.ndarray
        The wavelength grid in microns.
    """
    radmc_dir = Path(radmc_dir)
    wavelengths = np.logspace(np.log10(wave_min), np.log10(wave_max), n_wave)

    filepath = radmc_dir / "wavelength_micron.inp"
    with open(filepath, "w") as f:
        f.write(f"{len(wavelengths)}\n")
        for w in wavelengths:
            f.write(f"{w:e}\n")

    print(f"    wavelength_micron.inp ({n_wave} points, {wave_min}-{wave_max} um)")
    return wavelengths


# ============================================================
#  Write stars.inp (from derived stellar properties)
# ============================================================

def write_stars(radmc_dir, stars, wavelengths):
    """
    Write the RADMC-3D stars.inp file.

    Uses the derived stellar properties (position in metres, radius in
    solar radii, temperature in K, mass in solar masses) and converts
    to CGS as required by RADMC-3D.

    Parameters
    ----------
    radmc_dir : str or Path
        Output directory.
    stars : list of dict
        Stellar properties as returned by polaris.opacity.derive_stars_properties().
        Each dict has: pos_m, radius_rsun, mass_msun, temperature_k.
    wavelengths : np.ndarray
        Wavelength grid in microns.
    """
    radmc_dir = Path(radmc_dir)
    filepath = radmc_dir / "stars.inp"

    with open(filepath, "w") as f:
        f.write("2\n")
        f.write(f"{len(stars)} {len(wavelengths)}\n")

        for star in stars:
            r_cm = star["radius_rsun"] * R_sun        # R_sun is in cm
            m_g = star["mass_msun"] * M_sun            # M_sun is in g
            pos_cm = [p * 1e2 for p in star["pos_m"]]  # m -> cm

            f.write(
                f"{r_cm: .6e} {m_g: .6e} "
                f"{pos_cm[0]: .6e} {pos_cm[1]: .6e} {pos_cm[2]: .6e}\n"
            )

        for w in wavelengths:
            f.write(f"{w:e}\n")

        for star in stars:
            f.write(f"-{star['temperature_k']:e}\n")

    print(f"    stars.inp ({len(stars)} stars)")


# ============================================================
#  Write radmc3d.inp control file
# ============================================================

def write_radmc3d_control(radmc_dir, nphot=1_000_000, setthreads=8,
                          scattering_mode=1, modified_random_walk=1,
                          rto_style=3, rto_single=1):
    """
    Write the RADMC-3D radmc3d.inp control file.

    Parameters
    ----------
    radmc_dir : str or Path
        Output directory.
    nphot : int
        Number of photon packages for mctherm and scattering.
    setthreads : int
        Number of OpenMP threads.
    scattering_mode : int
        Scattering mode (1 = isotropic).
    modified_random_walk : int
        Enable modified random walk (1 = yes).
    rto_style : int
        Output style for dust_temperature (3 = binary).
    rto_single : int
        Single precision output (1 = yes).
    """
    radmc_dir = Path(radmc_dir)
    filepath = radmc_dir / "radmc3d.inp"

    with open(filepath, "w") as f:
        f.write(f"nphot = {int(nphot)}\n")
        f.write(f"nphot_scat = {int(nphot)}\n")
        f.write(f"setthreads = {setthreads}\n")
        f.write(f"scattering_mode = {scattering_mode}\n")
        f.write(f"scattering_mode_max = {scattering_mode}\n")
        f.write(f"modified_random_walk = {modified_random_walk}\n")
        f.write(f"rto_style = {rto_style}\n")
        f.write(f"rto_single = {rto_single}\n")

    print(f"    radmc3d.inp (nphot={int(nphot)}, threads={setthreads})")


# ============================================================
#  High-level: prepare all RADMC-3D inputs
# ============================================================

def prepare_radmc3d_inputs(
    ramses_dir,
    radmc_dir,
    polaris_data_dir,
    n_dust=None,
    f_acc=0.1,
    wave_min=0.27,
    wave_max=3000,
    n_wavelengths=200,
    nphot=1_000_000,
    setthreads=8,
    scattering_mode=1,
):
    """
    Full Step 5: convert POLARIS opacities and write all remaining
    RADMC-3D input files.

    This writes into radmc_dir:
        dustkappa_dust_mixture_001.inp, ..., dustkappa_dust_mixture_N.inp
        dustopac.inp
        wavelength_micron.inp
        stars.inp
        radmc3d.inp

    Note: amr_grid.inp and dust_density.inp are already written by
    Convert.to_radmc() (Step 3) and should already exist in radmc_dir.

    Parameters
    ----------
    ramses_dir : str or Path
        RAMSES output directory (used to derive stellar properties).
    radmc_dir : str or Path
        RADMC-3D working directory.
    polaris_data_dir : str or Path
        Path to the POLARIS data/ directory containing dust_mixture_*.dat.
    n_dust : int or None
        Number of dust species. If None, auto-detected from the number
        of dust_mixture_*.dat files in polaris_data_dir.
    f_acc : float
        Accretion efficiency for deriving stellar radii.
    wave_min, wave_max : float
        Wavelength range in microns.
    n_wavelengths : int
        Number of wavelength points.
    nphot : int
        Number of photon packages.
    setthreads : int
        Number of OpenMP threads.
    scattering_mode : int
        Scattering mode.

    Returns
    -------
    radmc_dir : Path
        The RADMC-3D directory (for chaining).
    """
    ramses_dir = Path(ramses_dir)
    radmc_dir = Path(radmc_dir)
    polaris_data_dir = Path(polaris_data_dir)

    radmc_dir.mkdir(parents=True, exist_ok=True)

    # Auto-detect n_dust from POLARIS output if not provided
    if n_dust is None:
        dust_files = sorted(polaris_data_dir.glob("dust_mixture_*.dat"))
        n_dust = len(dust_files)
        if n_dust == 0:
            raise FileNotFoundError(
                f"No dust_mixture_*.dat files found in {polaris_data_dir}. "
                "Run the POLARIS opacity step first (Step 4)."
            )
        print(f"Auto-detected {n_dust} dust species from POLARIS output.")

    print(f"\nPreparing RADMC-3D input files in: {radmc_dir}")
    print(f"POLARIS opacity data from: {polaris_data_dir}\n")

    # 1) Convert POLARIS opacities -> dustkappa_*.inp
    kappa_names = convert_polaris_opacities(
        polaris_data_dir, radmc_dir, n_dust,
    )

    # 2) dustopac.inp
    write_dustopac(radmc_dir, kappa_names)

    # 3) wavelength_micron.inp
    wavelengths = write_wavelength_micron(
        radmc_dir, wave_min=wave_min, wave_max=wave_max, n_wave=n_wavelengths,
    )

    # 4) stars.inp
    stars, boxlen_pc = derive_stars_properties(ramses_dir, f_acc=f_acc)
    if not stars:
        print("WARNING: No stars with luminosity data found. "
              "stars.inp will not be written.")
    else:
        write_stars(radmc_dir, stars, wavelengths)

    # 5) radmc3d.inp
    write_radmc3d_control(
        radmc_dir, nphot=nphot, setthreads=setthreads,
        scattering_mode=scattering_mode,
    )

    print(f"\nAll RADMC-3D input files prepared successfully in: {radmc_dir}")
    return radmc_dir
