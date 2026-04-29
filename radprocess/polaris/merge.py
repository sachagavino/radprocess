"""
_____________________________________________________________________________________________________________
file name: merge
last update: Apr 2026
language: > PYTHON 3.9
short description: Merge RADMC-3D dust temperatures into the POLARIS binary octree grid.

                   After POLARIS runs CMD_TEMP (Step 4), it produces a grid_temp.dat file
                   that contains the POLARIS grid with POLARIS-computed temperatures.
                   After RADMC-3D runs mctherm (Step 6), it produces dust_temperature.bdat.

                   This module replaces the POLARIS dust temperatures (parameter ID = 2)
                   in the POLARIS grid with the RADMC-3D values, producing a new grid file
                   (grid_temp.radmc3d.dat) that is physically consistent with the RADMC-3D
                   thermal solution. This merged grid is then used for the final POLARIS
                   imaging run (Step 8).
_____________________________________________________________________________________________________________
"""

import os
import struct
from pathlib import Path

import numpy as np

from radprocess.polaris.opacity import _find_info_file, _read_ndust


# ============================================================
#  Read RADMC-3D binary temperature file
# ============================================================

def read_radmc3d_temperature(temp_file, n_dust):
    """
    Read dust_temperature.bdat produced by RADMC-3D mctherm.

    The binary format has a 4-entry header (each 8 bytes = int64),
    followed by float32 temperature values. The values are stored
    as all cells for species 1, then all cells for species 2, etc.
    (Fortran column-major order).

    Parameters
    ----------
    temp_file : str or Path
        Path to dust_temperature.bdat.
    n_dust : int
        Number of dust species.

    Returns
    -------
    temperatures : np.ndarray, shape (n_cells, n_dust)
        Dust temperature per cell and per species, in K.
    """
    temp_file = Path(temp_file)

    with open(temp_file, "rb") as f:
        # Skip 4 header entries (each 8 bytes)
        f.seek(4 * 8)
        data = np.fromfile(f, dtype=np.float32)

    # Reshape: Fortran order means species vary slowest
    temperatures = data.reshape(-1, n_dust, order="F")

    return temperatures


# ============================================================
#  Merge temperatures into POLARIS grid
# ============================================================

def merge_temperature(
    polaris_grid_file,
    radmc_temp_file,
    output_file,
    n_dust,
    rename_original=True,
):
    """
    Replace dust temperature values (parameter ID = 2) in the POLARIS
    binary grid with temperatures from the RADMC-3D mctherm output.

    The POLARIS binary grid format is:
        Header:
            uint16  grid_id
            uint16  n_param
            uint16  param_id[0], ..., param_id[n_param-1]
            float64 grid_size
        Octree nodes (recursive):
            uint16  is_leaf
            uint16  level
            if is_leaf:
                float32 data[0], ..., data[n_param-1]
            else:
                8 child nodes (recursion)

    Parameters
    ----------
    polaris_grid_file : str or Path
        Path to the POLARIS grid file to read (grid_temp.dat or
        grid_temp.polaris.dat if already renamed).
    radmc_temp_file : str or Path
        Path to the RADMC-3D dust_temperature.bdat file.
    output_file : str or Path
        Path to write the merged grid file (grid_temp.radmc3d.dat).
    n_dust : int
        Number of dust species.
    rename_original : bool
        If True and polaris_grid_file is named grid_temp.dat, rename it
        to grid_temp.polaris.dat before processing (preserves the original).

    Returns
    -------
    output_file : Path
        Path to the written merged grid file.
    num_cells : int
        Number of leaf cells processed.
    """
    polaris_grid_file = Path(polaris_grid_file)
    radmc_temp_file = Path(radmc_temp_file)
    output_file = Path(output_file)

    # Handle renaming of original grid file
    if rename_original and polaris_grid_file.name == "grid_temp.dat":
        renamed = polaris_grid_file.parent / "grid_temp.polaris.dat"
        if renamed.exists():
            renamed.unlink()
        polaris_grid_file.rename(renamed)
        print(f"Renamed original grid to: {renamed}")
        polaris_grid_file = renamed

    if not polaris_grid_file.exists():
        raise FileNotFoundError(
            f"POLARIS grid file not found: {polaris_grid_file}"
        )

    # Load RADMC-3D temperatures
    print(f"Reading RADMC-3D temperatures from: {radmc_temp_file}")
    temperatures = read_radmc3d_temperature(radmc_temp_file, n_dust)
    num_leaf_cells = temperatures.shape[0]
    print(f"  Loaded temperatures for {num_leaf_cells} cells, {n_dust} dust species.")
    print(f"  T range: {temperatures.min():.1f} - {temperatures.max():.1f} K")

    # Process the POLARIS grid
    print(f"\nMerging into: {output_file}")

    leaf_index = 0

    with open(polaris_grid_file, "rb") as in_f, \
         open(output_file, "wb") as out_f:

        # --- Header ---
        # grid_id (uint16)
        grid_id = in_f.read(2)
        out_f.write(grid_id)

        # n_param (uint16)
        n_param_buf = in_f.read(2)
        out_f.write(n_param_buf)
        n_param = struct.unpack("H", n_param_buf)[0]

        # parameter IDs (uint16 * n_param)
        param_ids_buf = in_f.read(2 * n_param)
        param_ids = list(struct.unpack("H" * n_param, param_ids_buf))
        out_f.write(param_ids_buf)

        # Find temperature parameter indices (ID = 2 in POLARIS convention)
        # ID 2 = dust temperature
        temp_indices = [i for i, pid in enumerate(param_ids) if pid == 2]

        if len(temp_indices) == 0:
            raise ValueError(
                f"No temperature parameters (ID=2) found in POLARIS grid. "
                f"Parameter IDs are: {param_ids}. "
                f"The POLARIS opacity run may not have produced temperature data."
            )

        if len(temp_indices) != n_dust:
            raise ValueError(
                f"Mismatch: found {len(temp_indices)} temperature parameters (ID=2) "
                f"but expected {n_dust} dust species."
            )

        print(f"  Found {len(temp_indices)} temperature parameters at indices: {temp_indices}")

        # grid_size (float64)
        grid_size = in_f.read(8)
        out_f.write(grid_size)

        # --- Octree traversal ---
        while leaf_index < num_leaf_cells:
            # is_leaf (uint16)
            buf = in_f.read(2)
            if not buf:
                break
            out_f.write(buf)
            is_leaf = struct.unpack("H", buf)[0]

            # level (uint16)
            buf = in_f.read(2)
            if not buf:
                break
            out_f.write(buf)

            if is_leaf == 1:
                dust_temp_counter = 0
                for j in range(n_param):
                    val_buf = in_f.read(4)

                    if j in temp_indices:
                        # Replace with RADMC-3D temperature
                        temp_val = temperatures[leaf_index, dust_temp_counter]
                        out_f.write(struct.pack("f", temp_val))
                        dust_temp_counter += 1
                    else:
                        # Copy original value unchanged
                        out_f.write(val_buf)

                leaf_index += 1

            if leaf_index % 50000 == 0:
                print(f"  Processed {leaf_index}/{num_leaf_cells} cells...")

    print(f"  Processed {leaf_index} leaf cells total.")
    print(f"\nMerged grid written successfully: {output_file}\n")

    return output_file, leaf_index


# ============================================================
#  High-level orchestrator
# ============================================================

def merge_radmc3d_temperature(
    polaris_dir,
    radmc_dir,
    n_dust=None,
    ramses_dir=None,
):
    """
    Full Step 7: merge RADMC-3D dust temperatures into the POLARIS grid.

    Reads grid_temp.dat from the POLARIS directory (produced by Step 4),
    reads dust_temperature.bdat from the RADMC-3D directory (produced by
    Step 6), and writes grid_temp.radmc3d.dat back to the POLARIS directory.

    Parameters
    ----------
    polaris_dir : str or Path
        POLARIS working directory containing grid_temp.dat.
    radmc_dir : str or Path
        RADMC-3D working directory containing dust_temperature.bdat.
    n_dust : int or None
        Number of dust species. If None and ramses_dir is provided,
        auto-detected from the RAMSES info file.
    ramses_dir : str or Path or None
        RAMSES output directory (used to auto-detect n_dust).

    Returns
    -------
    merged_grid : Path
        Path to the merged grid file (grid_temp.radmc3d.dat).
    """
    polaris_dir = Path(polaris_dir)
    radmc_dir = Path(radmc_dir)

    # Auto-detect n_dust
    if n_dust is None:
        if ramses_dir is None:
            raise ValueError(
                "Either n_dust or ramses_dir must be provided "
                "so the number of dust species can be determined."
            )
        info_path = _find_info_file(ramses_dir)
        n_dust = _read_ndust(info_path)
        if n_dust == 0:
            n_dust = 1
            print("No dust in simulation. Assuming 1 virtual dust species.")
        print(f"Auto-detected {n_dust} dust species from RAMSES info.\n")

    # Locate input files
    # The POLARIS opacity run (Step 4) produces grid_temp.dat
    grid_temp = polaris_dir / "grid_temp.dat"
    grid_temp_renamed = polaris_dir / "grid_temp.polaris.dat"

    if grid_temp.exists():
        polaris_grid_file = grid_temp
    elif grid_temp_renamed.exists():
        polaris_grid_file = grid_temp_renamed
    else:
        raise FileNotFoundError(
            f"Neither grid_temp.dat nor grid_temp.polaris.dat found in "
            f"{polaris_dir}. Run the POLARIS opacity step first (Step 4)."
        )

    # RADMC-3D temperature
    radmc_temp_file = radmc_dir / "dust_temperature.bdat"
    if not radmc_temp_file.exists():
        radmc_temp_file = radmc_dir / "dust_temperature.dat"
    if not radmc_temp_file.exists():
        raise FileNotFoundError(
            f"No dust_temperature file found in {radmc_dir}. "
            "Run radmc3d mctherm first (Step 6)."
        )

    # Output
    output_file = polaris_dir / "grid_temp.radmc3d.dat"

    merged_grid, n_cells = merge_temperature(
        polaris_grid_file=polaris_grid_file,
        radmc_temp_file=radmc_temp_file,
        output_file=output_file,
        n_dust=n_dust,
    )

    return merged_grid