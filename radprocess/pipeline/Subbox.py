"""
_____________________________________________________________________________________________________________
file name: Subbox
last update: June 2026
language: > PYTHON 3.9
short description: Extract subboxes around individual sinks from a RAMSES Zarr grid
                   and produce one RADMC-3D or POLARIS folder per isolated sink.
_____________________________________________________________________________________________________________
"""

import sys
import struct
from pathlib import Path

import numpy as np

from radprocess.pipeline.OcTree import OcTree, CellOct
from radprocess import radmc3d
from radprocess import ramses
from radprocess.constants.constants import (
    Ggram, kbol, pc2cm, pc2m, amu, au2m, au2cm,
    M_sun, L_sun, R_sun, sigma,
)


# ============================================================
#  Sink filtering
# ============================================================

def filter_sinks(sinks, isolation_radius_au=100.0, require_luminosity=True):
    """
    From a SinkInfo object, return indices of isolated sinks suitable
    for individual subbox post-processing.

    Algorithm
    ---------
    1. Optionally discard sinks that lack valid Lint and Teff
       (i.e. acc_lum <= 0 AND int_lum <= 0).
    2. Group remaining sinks by proximity: two sinks closer than
       `isolation_radius_au` belong to the same group.
       Within each group, keep only the most massive sink.

    Parameters
    ----------
    sinks : SinkInfo
        As returned by ramses.read.sink_info().
    isolation_radius_au : float
        Minimum separation (in AU) for a sink to be considered isolated.
        For groups of sinks closer than this, only the most massive is kept.
    require_luminosity : bool
        If True, discard sinks with acc_lum <= 0 AND int_lum <= 0.

    Returns
    -------
    keep_indices : np.ndarray of int
        Row indices (into sinks.data) of the sinks to keep.
    """
    cols = sinks.columns
    x_idx = cols.index("x")
    y_idx = cols.index("y")
    z_idx = cols.index("z")
    m_idx = cols.index("M[Msol]")
    acclum_idx = cols.index("acc_lum[Lsol]")
    intlum_idx = cols.index("int_lum[Lsol]")

    data = sinks.data
    n = data.shape[0]

    # ----------------------------------------------------------
    # Step 1: luminosity filter
    # ----------------------------------------------------------
    if require_luminosity:
        has_lum = (data[:, acclum_idx] > 0) | (data[:, intlum_idx] > 0)
        candidates = np.where(has_lum)[0]
    else:
        candidates = np.arange(n)

    if len(candidates) == 0:
        return np.array([], dtype=int)

    # ----------------------------------------------------------
    # Step 2: greedy grouping by proximity (friends-of-friends)
    # ----------------------------------------------------------
    pos_pc = data[candidates][:, [x_idx, y_idx, z_idx]]
    pos_au = pos_pc * (pc2m / au2m)  # pc -> AU

    masses = data[candidates, m_idx]

    visited = np.zeros(len(candidates), dtype=bool)
    keep = []

    for i in range(len(candidates)):
        if visited[i]:
            continue

        dx = pos_au[:, 0] - pos_au[i, 0]
        dy = pos_au[:, 1] - pos_au[i, 1]
        dz = pos_au[:, 2] - pos_au[i, 2]
        dist = np.sqrt(dx**2 + dy**2 + dz**2)

        group = np.where((dist < isolation_radius_au) & (~visited))[0]
        visited[group] = True

        best = group[np.argmax(masses[group])]
        keep.append(candidates[best])

    return np.array(keep, dtype=int)


# ============================================================
#  Shared subbox extraction
# ============================================================

def _extract_subbox(root, sink_idx, sinks, box_half_width_au, hole_au,
                    f_acc, boxlen):
    """
    Extract cells from the Zarr within a box around a sink, recenter,
    relevel, and compute sink stellar properties.

    Parameters
    ----------
    boxlen : float
        RAMSES boxlen from the info file (dimensionless code units).
        Sink positions in the CSV are in [0, boxlen].

    Returns a dict with all the data needed by format-specific writers.
    """
    l_m = root.attrs.get("l_m")
    l_cm = root.attrs.get("l_cm")

    level = np.array(root["level"])
    x = np.array(root["x"])
    y = np.array(root["y"])
    z = np.array(root["z"])

    # Dust density (always needed)
    if "dust_massdensities" in root:
        dust_massdensity = np.array(root["dust_massdensities"])
    elif "dust_massdensity" in root:
        dust_massdensity = np.array(root["dust_massdensity"])
        if dust_massdensity.ndim == 1:
            dust_massdensity = dust_massdensity[:, np.newaxis]
    else:
        raise RuntimeError("No dust density field found in Zarr.")

    # Derive nb_species from actual array shape, not from attrs
    nb_species = dust_massdensity.shape[1]

    # Gas density (needed for POLARIS)
    gas_massdensity = None
    if "gas_massdensity" in root:
        gas_massdensity = np.array(root["gas_massdensity"])

    # Temperature (needed for POLARIS)
    Tgas = None
    if "Tgas" in root:
        Tgas = np.array(root["Tgas"])
    elif "Tdust" in root:
        Tgas = np.array(root["Tdust"])

    # B-field (needed for POLARIS)
    Bx = By = Bz = None
    if "Bx" in root:
        Bx = np.array(root["Bx"])
        By = np.array(root["By"])
        Bz = np.array(root["Bz"])

    # Velocity field (needed for POLARIS line RT)
    Vel = None
    if "velocity" in root:
        Vel = np.array(root["velocity"])  # shape (N, 3), in m/s

    # ----------------------------------------------------------
    # Sink properties
    # ----------------------------------------------------------
    cols = sinks.columns
    x_col = cols.index("x")
    y_col = cols.index("y")
    z_col = cols.index("z")
    m_col = cols.index("M[Msol]")
    accrate_col = cols.index("acc_rate[Msol/y]")
    acclum_col = cols.index("acc_lum[Lsol]")
    intlum_col = cols.index("int_lum[Lsol]")

    sec_yr = 365.25 * 24 * 3600

    sink_pos_raw = sinks.data[sink_idx, [x_col, y_col, z_col]]
    # Sink positions are in RAMSES code units [0, boxlen].
    # Normalize to [0, 1], then scale to physical [0, l_m] and center.
    # This matches how cell positions are stored: cells.points * l_m.
    sink_pos_m = (sink_pos_raw / boxlen - 0.5) * l_m

    sink_mass_g = sinks.data[sink_idx, m_col] * M_sun
    acc_rate_g = sinks.data[sink_idx, accrate_col] * M_sun / sec_yr
    lint = sinks.data[sink_idx, intlum_col] * L_sun
    lacc = sinks.data[sink_idx, acclum_col] * L_sun
    ltot = lint + lacc

    if lacc > 0:
        sink_radius = f_acc * Ggram * sink_mass_g * acc_rate_g / lacc
        sink_teff = (ltot / (4 * np.pi * sigma * sink_radius**2))**0.25
    else:
        sink_radius = 1.0 * R_sun
        sink_teff = 5e3

    # ----------------------------------------------------------
    # Compute box size and snap center to AMR grid
    # ----------------------------------------------------------
    # Architecture: the grid layer extracts a PADDED, GRID-ALIGNED box
    # (artifact-free octree), while the imaging layer (RADMC-3D regrid,
    # POLARIS detector) crops to the user's requested FOV centered on
    # the sink.
    #
    # The box uses k_extract = k_base - 1 (one level coarser = 2× bigger).
    # The origin snaps at level k_base (finer grid), giving a maximum
    # sink offset of l_snap/2 per axis.  Since l_box = 2·l_snap and
    # l_snap ≥ 2·hw (by construction), the margin from the sink to the
    # box edge is always ≥ l_snap/2 ≥ hw.
    # ----------------------------------------------------------
    hw_m = box_half_width_au * au2m

    cx = x - 0.5 * l_m
    cy = y - 0.5 * l_m
    cz = z - 0.5 * l_m

    # k_base: the smallest k such that l_m/2^k >= 2*hw
    l_box_raw = 2.0 * hw_m
    k_base = int(np.floor(np.log2(l_m / l_box_raw)))
    if k_base < 0:
        k_base = 0
    l_snap = l_m / (2**k_base)   # snap resolution

    # k_extract: one level coarser (2× bigger) for padding
    k_extract = max(k_base - 1, 0)
    l_box = l_m / (2**k_extract)
    hw_actual = l_box / 2.0
    level_offset = k_extract
    half_domain = l_m / 2.0

    # Snap the box origin at level k_base (finer than the box size).
    # This ensures all cells at level ≥ k_base fall on exact octree
    # subdivision points, while giving finer positioning than snapping
    # at k_extract.
    desired_origin_x = sink_pos_m[0] - hw_actual
    desired_origin_y = sink_pos_m[1] - hw_actual
    desired_origin_z = sink_pos_m[2] - hw_actual

    n_x = int(np.round((desired_origin_x + half_domain) / l_snap))
    n_y = int(np.round((desired_origin_y + half_domain) / l_snap))
    n_z = int(np.round((desired_origin_z + half_domain) / l_snap))

    # The box spans n_span snap cells
    n_span = int(round(l_box / l_snap))  # = 2
    n_snap_max = int(round(l_m / l_snap)) - n_span
    n_x = max(0, min(n_x, n_snap_max))
    n_y = max(0, min(n_y, n_snap_max))
    n_z = max(0, min(n_z, n_snap_max))

    box_origin_x = -half_domain + n_x * l_snap
    box_origin_y = -half_domain + n_y * l_snap
    box_origin_z = -half_domain + n_z * l_snap

    snap_center_x = box_origin_x + hw_actual
    snap_center_y = box_origin_y + hw_actual
    snap_center_z = box_origin_z + hw_actual

    # ----------------------------------------------------------
    # Extract cells that fall within the box
    # ----------------------------------------------------------
    mask = (
        (cx >= box_origin_x) & (cx < box_origin_x + l_box) &
        (cy >= box_origin_y) & (cy < box_origin_y + l_box) &
        (cz >= box_origin_z) & (cz < box_origin_z + l_box)
    )

    nb_cells_sub = int(mask.sum())
    if nb_cells_sub == 0:
        raise RuntimeError(
            f"No cells found in subbox of +/-{box_half_width_au} AU "
            f"around sink {sink_idx}. Check coordinates."
        )

    sub_cx = cx[mask]
    sub_cy = cy[mask]
    sub_cz = cz[mask]
    sub_level = level[mask]
    sub_dust = dust_massdensity[mask]
    sub_gas = gas_massdensity[mask] if gas_massdensity is not None else None
    sub_T = Tgas[mask] if Tgas is not None else None
    sub_Bx = Bx[mask] if Bx is not None else None
    sub_By = By[mask] if By is not None else None
    sub_Bz = Bz[mask] if Bz is not None else None
    sub_Vel = Vel[mask] if Vel is not None else None

    # ----------------------------------------------------------
    # Re-level: if coarse cells have level < k_extract, enlarge the
    # box and re-snap (grid-aligned at the coarser level).
    # ----------------------------------------------------------
    sub_level_local = (sub_level - level_offset).astype(int)

    if sub_level_local.min() < 0:
        min_level_full = int(sub_level.min())
        level_offset = min_level_full
        k_extract = min_level_full
        l_box = l_m / (2**k_extract)
        hw_actual = l_box / 2.0
        sub_level_local = (sub_level - level_offset).astype(int)

        # Re-snap with the larger box (grid-aligned at k_extract)
        desired_origin_x = sink_pos_m[0] - hw_actual
        desired_origin_y = sink_pos_m[1] - hw_actual
        desired_origin_z = sink_pos_m[2] - hw_actual

        n_x = int(np.round((desired_origin_x + half_domain) / l_box))
        n_y = int(np.round((desired_origin_y + half_domain) / l_box))
        n_z = int(np.round((desired_origin_z + half_domain) / l_box))

        n_max = 2**k_extract - 1
        n_x = max(0, min(n_x, n_max))
        n_y = max(0, min(n_y, n_max))
        n_z = max(0, min(n_z, n_max))

        box_origin_x = -half_domain + n_x * l_box
        box_origin_y = -half_domain + n_y * l_box
        box_origin_z = -half_domain + n_z * l_box

        snap_center_x = box_origin_x + hw_actual
        snap_center_y = box_origin_y + hw_actual
        snap_center_z = box_origin_z + hw_actual

        mask = (
            (cx >= box_origin_x) & (cx < box_origin_x + l_box) &
            (cy >= box_origin_y) & (cy < box_origin_y + l_box) &
            (cz >= box_origin_z) & (cz < box_origin_z + l_box)
        )
        sub_cx = cx[mask]
        sub_cy = cy[mask]
        sub_cz = cz[mask]
        sub_level = level[mask]
        sub_level_local = (sub_level - level_offset).astype(int)
        sub_dust = dust_massdensity[mask]
        sub_gas = gas_massdensity[mask] if gas_massdensity is not None else None
        sub_T = Tgas[mask] if Tgas is not None else None
        sub_Bx = Bx[mask] if Bx is not None else None
        sub_By = By[mask] if By is not None else None
        sub_Bz = Bz[mask] if Bz is not None else None
        sub_Vel = Vel[mask] if Vel is not None else None
        nb_cells_sub = int(mask.sum())

    # Cell positions relative to box origin (for octree insertion)
    rel_x = sub_cx - box_origin_x
    rel_y = sub_cy - box_origin_y
    rel_z = sub_cz - box_origin_z

    # Sink offset from box center (non-zero due to grid snapping)
    sink_offset_x = sink_pos_m[0] - snap_center_x
    sink_offset_y = sink_pos_m[1] - snap_center_y
    sink_offset_z = sink_pos_m[2] - snap_center_z

    max_level_local = int(sub_level_local.max())
    min_level_local = int(sub_level_local.min())

    offset_au = np.sqrt(sink_offset_x**2 + sink_offset_y**2 + sink_offset_z**2) / au2m
    margin_au = (hw_actual - max(abs(sink_offset_x), abs(sink_offset_y),
                                  abs(sink_offset_z))) / au2m - box_half_width_au

    print(f"\n--- Subbox for sink {sink_idx} ---")
    print(f"    Sink position (AU): "
          f"({sink_pos_m[0]/au2m:.1f}, {sink_pos_m[1]/au2m:.1f}, {sink_pos_m[2]/au2m:.1f})")
    print(f"    Box center snapped (AU): "
          f"({snap_center_x/au2m:.1f}, {snap_center_y/au2m:.1f}, {snap_center_z/au2m:.1f})")
    print(f"    Sink offset from center (AU): "
          f"({sink_offset_x/au2m:.1f}, {sink_offset_y/au2m:.1f}, {sink_offset_z/au2m:.1f})")
    print(f"    Requested half-width: {box_half_width_au} AU")
    print(f"    Actual subbox side:   {l_box/au2m:.1f} AU "
          f"(padded for grid alignment)")
    print(f"    Margin to box edge:   {margin_au:.1f} AU")
    print(f"    Cells in subbox:      {nb_cells_sub}")
    print(f"    Level range (local):  {min_level_local} - {max_level_local}")
    print(f"    Level offset:         {level_offset}")
    print(f"    Dust species:         {nb_species}\n")

    return {
        "l_m": l_m,
        "l_box": l_box,
        "nb_cells_sub": nb_cells_sub,
        "nb_species": nb_species,
        "max_level_local": max_level_local,
        "min_level_local": min_level_local,
        # Recentered cell positions
        "rel_x": rel_x,
        "rel_y": rel_y,
        "rel_z": rel_z,
        # Original centered positions (for hole digging)
        "sub_cx": sub_cx,
        "sub_cy": sub_cy,
        "sub_cz": sub_cz,
        "sub_level_local": sub_level_local,
        # Physical data
        "sub_dust": sub_dust,
        "sub_gas": sub_gas,
        "sub_T": sub_T,
        "sub_Bx": sub_Bx,
        "sub_By": sub_By,
        "sub_Bz": sub_Bz,
        "sub_Vel": sub_Vel,
        # Sink properties
        "sink_pos_m": sink_pos_m,
        "sink_mass_g": sink_mass_g,
        "sink_radius": sink_radius,
        "sink_teff": sink_teff,
        "hole_au": hole_au,
        # Sink offset from box center (for star positioning and imaging)
        "sink_offset_m": np.array([sink_offset_x, sink_offset_y, sink_offset_z]),
        # User's requested half-width (for downstream cropping)
        "requested_hw_au": box_half_width_au,
    }


# ============================================================
#  RADMC-3D subbox writer
# ============================================================

def build_subbox_radmc(
    root, ramses_dir, output_dir, sink_idx, sinks,
    box_half_width_au, hole_au, f_acc, boxlen,
    gridstyle="octtree", coordsystem="cartesian",
    lmin=1e-3, lmax=1e4, nlam=210,
):
    """
    Extract a subbox around a sink and write RADMC-3D files.

    Parameters
    ----------
    root : zarr.Group
    ramses_dir : str or Path
    output_dir : str or Path
    sink_idx : int
    sinks : SinkInfo
    box_half_width_au : float
    hole_au : float
    f_acc : float
    boxlen : float
    gridstyle, coordsystem : str
    lmin, lmax : float
    nlam : int

    Returns
    -------
    grid : list
    densityarray : np.ndarray
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sub = _extract_subbox(root, sink_idx, sinks, box_half_width_au,
                          hole_au, f_acc, boxlen)

    l_box = sub["l_box"]
    l_box_cm = l_box * 100.0
    nb_cells_sub = sub["nb_cells_sub"]
    nb_species = sub["nb_species"]


    # Vectorized hole digging
    if hole_au > 0:
        hole_radius2 = (hole_au * au2m)**2
        d2 = ((sub["sub_cx"] - sub["sink_pos_m"][0])**2 +
              (sub["sub_cy"] - sub["sink_pos_m"][1])**2 +
              (sub["sub_cz"] - sub["sink_pos_m"][2])**2)
        hole_mask = d2 <= hole_radius2
        sub["sub_dust"][hole_mask, :] = 0.0
        print(f"    Hole digging: zeroed {hole_mask.sum()} cells")

    # Build octree
    tree = OcTree(0.0, 0.0, 0.0, l_box)
    tree.nr_of_cells = nb_cells_sub

    print("Constructing subbox octree (RADMC-3D)...")
    for i in range(nb_cells_sub):
        cell = CellOct(sub["rel_x"][i], sub["rel_y"][i], sub["rel_z"][i],
                        0, int(sub["sub_level_local"][i]))
        cell.data = sub["sub_dust"][i, :].tolist()
        tree.insertInTree(tree.root, cell, 0)

        if i % 10000 == 0:
            sys.stdout.write(f"    Constructing octree: {100.*i/nb_cells_sub:.1f}%\r")
            sys.stdout.flush()

    print("    Constructing octree: done\n")

    tree.reset_counter()
    check = tree.checkOcTree(tree.root)
    if not check:
        raise RuntimeError(f"Octree integrity check failed for sink {sink_idx}!")
    print("    Octree structure: OK\n")

    # Write RADMC-3D files
    tree.cell_counter = 0
    grid = []
    density = []
    tree._n_species = nb_species
    tree.writeOcTree_radmc(tree.root, grid, density)
    densityarray = np.array(density)

    # The actual number of leaves may be larger than nb_cells_sub
    # because empty branch nodes are padded with zero-density leaves
    nb_leaves_actual = len(densityarray)

    print(f"    Leaf cells (from Zarr): {nb_cells_sub}")
    print(f"    Leaf cells (with padding): {nb_leaves_actual}")

    print("    Writing amr_grid.inp ...")
    radmc3d.write.amr_grid(
        output_dir, grid, sub["max_level_local"], nb_leaves_actual, l_box_cm,
        gridstyle=gridstyle, coordsystem=coordsystem,
    )

    print("    Writing dust_density.inp ...")
    radmc3d.write.dust_density(
        output_dir, densityarray, nb_leaves_actual, nb_species,
        gridstyle=gridstyle,
    )

    print(f"    Sink {sink_idx} done (RADMC-3D) -> {output_dir}\n")
    print(f"    NOTE: stars.inp will be written by prepare_radmc3d_inputs(subbox=True)")

    # Save sink offset for stars.inp positioning
    offset_file = output_dir / "sink_offset.txt"
    offset_cm = sub["sink_offset_m"] * 100.0  # m -> cm (RADMC-3D uses CGS)
    with open(offset_file, "w") as f:
        f.write(f"{offset_cm[0]:.10e} {offset_cm[1]:.10e} {offset_cm[2]:.10e}\n")

    # Save the user's requested FOV half-width for RADMC-3D subbox_regrid
    fov_file = output_dir / "requested_hw_au.txt"
    with open(fov_file, "w") as f:
        f.write(f"{sub['requested_hw_au']:.6e}\n")

    return grid, densityarray


# ============================================================
#  POLARIS subbox writer
# ============================================================

def build_subbox_polaris(
    root, ramses_dir, output_dir, sink_idx, sinks,
    box_half_width_au, hole_au, f_acc, boxlen,
):
    """
    Extract a subbox around a sink and write a POLARIS binary grid file.

    Parameters
    ----------
    root : zarr.Group
    ramses_dir : str or Path
    output_dir : str or Path
    sink_idx : int
    sinks : SinkInfo
    box_half_width_au : float
    hole_au : float
    f_acc : float
    boxlen : float

    Returns
    -------
    output_file : Path
    """
    POLARIS_GRID_ID = 20  # octree

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sub = _extract_subbox(root, sink_idx, sinks, box_half_width_au,
                          hole_au, f_acc, boxlen)

    l_box = sub["l_box"]
    nb_cells_sub = sub["nb_cells_sub"]
    nb_species = sub["nb_species"]
    n_dust = sub["sub_dust"].shape[1]
    hole_radius2 = (hole_au * au2m)**2

    # Determine available fields
    has_gas = sub["sub_gas"] is not None
    has_T = sub["sub_T"] is not None
    has_B = sub["sub_Bx"] is not None
    has_V = sub["sub_Vel"] is not None

    if not has_gas:
        print("WARNING: No gas density in Zarr. Using sum of dust densities.")
    if not has_T:
        print("WARNING: No temperature in Zarr. Setting T=10 K everywhere.")
    if not has_B:
        print("WARNING: No B-field in Zarr. Setting B=0.")
    if not has_V:
        print("WARNING: No velocity field in Zarr. Setting V=0.")

    # Build octree with POLARIS cell data:
    # Bx, By, Bz, Vx, Vy, Vz, gas_density (kg/m^3), Tgas (K), dust_density_1, ..., dust_density_N (kg/m^3)
    tree = OcTree(0.0, 0.0, 0.0, l_box)
    tree.nr_of_cells = nb_cells_sub

    print("Constructing subbox octree (POLARIS)...")
    for i in range(nb_cells_sub):
        # Densities: Zarr stores g/cc, POLARIS expects kg/m^3
        cell_dust = sub["sub_dust"][i, :].copy() * 1e3  # g/cc -> kg/m^3
        cell_gas = (sub["sub_gas"][i] * 1e3) if has_gas else cell_dust.sum()
        cell_T = sub["sub_T"][i] if has_T else 10.0
        cell_Bx = float(sub["sub_Bx"][i]) if has_B else 0.0
        cell_By = float(sub["sub_By"][i]) if has_B else 0.0
        cell_Bz = float(sub["sub_Bz"][i]) if has_B else 0.0
        # Velocity: Zarr stores m/s, POLARIS expects m/s (SI)
        cell_Vx = float(sub["sub_Vel"][i, 0]) if has_V else 0.0
        cell_Vy = float(sub["sub_Vel"][i, 1]) if has_V else 0.0
        cell_Vz = float(sub["sub_Vel"][i, 2]) if has_V else 0.0

        # Dig hole
        dx = sub["sub_cx"][i] - sub["sink_pos_m"][0]
        dy = sub["sub_cy"][i] - sub["sink_pos_m"][1]
        dz = sub["sub_cz"][i] - sub["sink_pos_m"][2]
        if dx**2 + dy**2 + dz**2 <= hole_radius2:
            cell_gas = 0.0
            cell_dust[:] = 0.0

        cell_data = [
            cell_Bx, cell_By, cell_Bz,
            cell_Vx, cell_Vy, cell_Vz,
            float(cell_gas),
            float(cell_T),
        ] + cell_dust.tolist()

        cell = CellOct(sub["rel_x"][i], sub["rel_y"][i], sub["rel_z"][i],
                        0, int(sub["sub_level_local"][i]))
        cell.data = cell_data
        tree.insertInTree(tree.root, cell, 0)

        if i % 10000 == 0:
            sys.stdout.write(f"    Constructing octree: {100.*i/nb_cells_sub:.1f}%\r")
            sys.stdout.flush()

    print("    Constructing octree: done\n")

    tree.reset_counter()
    check = tree.checkOcTree(tree.root)
    if not check:
        raise RuntimeError(f"Octree integrity check failed for sink {sink_idx}!")
    print("    Octree structure: OK\n")

    # Write POLARIS binary grid
    output_file = output_dir / f"ramses_grid_sink_{sink_idx:04d}.dat"

    # Data IDs: 4=Bx, 5=By, 6=Bz, 7=Vx, 8=Vy, 9=Vz, 28=gas density, 3=gas temperature, 29=dust density
    data_ids = [4, 5, 6, 7, 8, 9, 28, 3] + [29] * n_dust
    data_len = len(data_ids)

    print(f"    Writing POLARIS binary grid to: {output_file}")

    with open(output_file, "wb") as f:
        f.write(struct.pack("H", POLARIS_GRID_ID))
        f.write(struct.pack("H", data_len))
        for d_id in data_ids:
            f.write(struct.pack("H", d_id))
        f.write(struct.pack("d", l_box))  # grid size in meters

        tree.cell_counter = 0
        tree._n_data = data_len  # for empty leaf padding in subboxes
        tree.writeOcTree(f, tree.root)

    print(f"    Sink {sink_idx} done (POLARIS) -> {output_file}\n")

    # Save sink offset for source positioning (in metres for POLARIS)
    offset_file = output_dir / "sink_offset.txt"
    offset_m = sub["sink_offset_m"]
    with open(offset_file, "w") as f:
        f.write(f"{offset_m[0]:.10e} {offset_m[1]:.10e} {offset_m[2]:.10e}\n")

    # Save the user's requested FOV half-width
    fov_file = output_dir / "requested_hw_au.txt"
    with open(fov_file, "w") as f:
        f.write(f"{sub['requested_hw_au']:.6e}\n")

    return output_file