"""
_____________________________________________________________________________________________________________
file name: run
last update: Apr 2026
language: > PYTHON 3.9
short description: Execute RADMC-3D commands (mctherm, mcmono, image) from a given
                   working directory, with real-time stdout streaming and log capture.
_____________________________________________________________________________________________________________
"""

import sys
import subprocess
from pathlib import Path


def mctherm(radmc_dir, log_path=None, radmc3d_binary="radmc3d"):
    """
    Run ``radmc3d mctherm`` in the given directory.

    All RADMC-3D input files (amr_grid.inp, dust_density.inp,
    dustopac.inp, dustkappa_*.inp, stars.inp, wavelength_micron.inp,
    radmc3d.inp) must already exist in radmc_dir.

    The output is ``dust_temperature.bdat`` (or ``.dat`` depending on
    rto_style) written into the same directory.

    Parameters
    ----------
    radmc_dir : str or Path
        Working directory containing all RADMC-3D input files.
    log_path : str or Path or None
        If provided, stdout/stderr is also written to this log file.
    radmc3d_binary : str
        Name or path of the RADMC-3D executable.

    Returns
    -------
    temp_file : Path
        Path to the dust_temperature output file.

    Raises
    ------
    FileNotFoundError
        If required input files are missing.
    subprocess.CalledProcessError
        If RADMC-3D exits with a nonzero return code.
    """
    radmc_dir = Path(radmc_dir)

    # Verify essential input files exist
    required = [
        "amr_grid.inp",
        "dust_density.inp",
        "dustopac.inp",
        "radmc3d.inp",
        "wavelength_micron.inp",
    ]
    missing = [f for f in required if not (radmc_dir / f).exists()]
    if missing:
        raise FileNotFoundError(
            f"Missing RADMC-3D input files in {radmc_dir}: {missing}. "
            "Run the previous pipeline steps first."
        )

    # Check that at least one dustkappa file exists
    kappa_files = list(radmc_dir.glob("dustkappa_*.inp"))
    if not kappa_files:
        raise FileNotFoundError(
            f"No dustkappa_*.inp files found in {radmc_dir}. "
            "Run prepare_radmc3d_inputs() first (Step 5)."
        )

    command = [radmc3d_binary, "mctherm"]

    print(f"\nRunning: {' '.join(command)}")
    print(f"Working directory: {radmc_dir}")
    if log_path:
        print(f"Log file: {log_path}\n")

    process = subprocess.Popen(
        command,
        cwd=str(radmc_dir),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
    )

    log_file = open(log_path, "w") if log_path else None
    output_lines = []

    try:
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            output_lines.append(line)
            if log_file:
                log_file.write(line.rstrip("\r\n") + "\n")
    finally:
        if log_file:
            log_file.close()

    process.wait()

    # Check for errors (RADMC-3D may exit with code 0 even on abort)
    output_text = "".join(output_lines)
    if process.returncode != 0 or "ERROR" in output_text or "ABORTED" in output_text:
        raise RuntimeError(
            f"RADMC-3D mctherm failed in {radmc_dir}. "
            f"Check the log file for details."
        )

    # Find the output temperature file
    temp_bdat = radmc_dir / "dust_temperature.bdat"
    temp_dat = radmc_dir / "dust_temperature.dat"

    if temp_bdat.exists():
        temp_file = temp_bdat
    elif temp_dat.exists():
        temp_file = temp_dat
    else:
        raise RuntimeError(
            f"RADMC-3D completed but no dust_temperature file found in {radmc_dir}."
        )

    print(f"\nRADMC-3D mctherm completed successfully.")
    print(f"Temperature file: {temp_file}\n")

    return temp_file