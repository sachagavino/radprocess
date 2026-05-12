"""
_____________________________________________________________________________________________________________
file name: Subbox
last update: Apr 2026
language: > PYTHON 3.9
short description: Extract subboxes around individual sinks from a RAMSES Zarr grid
                   and produce one RADMC-3D folder per isolated sink.
_____________________________________________________________________________________________________________
"""

import sys
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
    # Positions in AU for the distance comparison
    pos_pc = data[candidates][:, [x_idx, y_idx, z_idx]]
    pos_au = pos_pc * (pc2m / au2m)  # pc -> AU

    masses = data[candidates, m_idx]

    visited = np.zeros(len(candidates), dtype=bool)
    keep = []

    for i in range(len(candidates)):
        if visited[i]:
            continue

        # Find all neighbors within isolation_radius
        dx = pos_au[:, 0] - pos_au[i, 0]
        dy = pos_au[:, 1] - pos_au[i, 1]
        dz = pos_au[:, 2] - pos_au[i, 2]
        dist = np.sqrt(dx**2 + dy**2 + dz**2)

        group = np.where((dist < isolation_radius_au) & (~visited))[0]
        visited[group] = True

        # Keep the most massive in the group
        best = group[np.argmax(masses[group])]
        keep.append(candidates[best])

    return np.array(keep, dtype=int)


# ============================================================
#  Subbox octree builder
# ============================================================

def build_subbox_radmc(
    root,
    ramses_dir,
    output_dir,
    sink_idx,
    sinks,
    box_half_width_au,
    hole_au,
    f_acc,
    boxlen_pc,
    gridstyle="octtree",
    coordsystem="cartesian",
    lmin=1e-3,
    lmax=1e4,
    nlam=210,
):
    """
    Extract a cubic subbox of side 2*box_half_width_au (in AU) centered
    on sink number `sink_idx`, build an octree from the cells that fall
    inside, and write a complete RADMC-3D folder.

    Parameters
    ----------
    root : zarr.Group
        The full-cloud Zarr grid (as produced by Convert.ramses()).
    ramses_dir : str or Path
        Path to the RAMSES output (needed only for context, not re-read).
    output_dir : str or Path
        Directory where RADMC-3D files will be written.
    sink_idx : int
        Row index of the target sink in sinks.data.
    sinks : SinkInfo
        Full sink table.
    box_half_width_au : float
        Half-width of the subbox in AU.
    hole_au : float
        Hole radius around the central sink in AU (density set to 0).
    f_acc : float
        Accretion efficiency factor for computing stellar radii.
    boxlen_pc : float
        RAMSES simulation box length in parsecs.
    gridstyle : str
        "octtree" (default) or "regular".
    coordsystem : str
        "cartesian" (default).
    lmin, lmax : float
        Wavelength range in microns for stars.inp.
    nlam : int
        Number of wavelength points.

    Returns
    -------
    grid : list
        The RADMC-3D octree grid array.
    densityarray : np.ndarray
        The RADMC-3D density array.
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ----------------------------------------------------------
    # 0) Read Zarr arrays
    # ----------------------------------------------------------
    l_m = root.attrs.get("l_m")
    l_cm = root.attrs.get("l_cm")
    nb_cells_total = root.attrs.get("nb_cells")
    nb_species = root.attrs.get("nb_species", 1)

    level = np.array(root["level"])
    x = np.array(root["x"])
    y = np.array(root["y"])
    z = np.array(root["z"])

    if "dust_massdensities" in root:
        dust_massdensity = np.array(root["dust_massdensities"])
    elif "dust_massdensity" in root:
        dust_massdensity = np.array(root["dust_massdensity"])[:, np.newaxis]

    # ----------------------------------------------------------
    # 1) Sink properties
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

    # Sink position in RAMSES code units (pc), then centered and in meters
    sink_pos_pc = sinks.data[sink_idx, [x_col, y_col, z_col]]
    sink_pos_centered_pc = sink_pos_pc - boxlen_pc / 2.0
    sink_pos_m = sink_pos_centered_pc * pc2m  # meters, centered

    sink_mass_g = sinks.data[sink_idx, m_col] * M_sun
    acc_rate_g = sinks.data[sink_idx, accrate_col] * M_sun / sec_yr
    lint = sinks.data[sink_idx, intlum_col] * L_sun
    lacc = sinks.data[sink_idx, acclum_col] * L_sun
    ltot = lint + lacc

    if lacc > 0:
        sink_radius = f_acc * Ggram * sink_mass_g * acc_rate_g / lacc
        sink_teff = (ltot / (4 * np.pi * sigma * sink_radius**2))**0.25
    else:
        sink_radius = 1.0 * R_sun  # fallback
        sink_teff = 5e3

    # ----------------------------------------------------------
    # 2) Spatial mask: select cells inside the subbox
    # ----------------------------------------------------------
    hw_m = box_half_width_au * au2m

    # Cell positions in the Zarr are in meters from the domain corner (0 to l_m).
    # Centered positions:
    cx = x - 0.5 * l_m
    cy = y - 0.5 * l_m
    cz = z - 0.5 * l_m

    mask = (
        (np.abs(cx - sink_pos_m[0]) <= hw_m) &
        (np.abs(cy - sink_pos_m[1]) <= hw_m) &
        (np.abs(cz - sink_pos_m[2]) <= hw_m)
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

    # ----------------------------------------------------------
    # 3) Re-level: the octree base cell must be level 0 relative
    #    to the subbox domain size.
    #
    #    In the full domain: cell size at level L = l_m / 2^L
    #    In the subbox:      cell size at level L' = l_box / 2^L'
    #    where l_box = 2 * hw_m.
    #
    #    Matching: l_m / 2^L = l_box / 2^L'
    #    => L' = L - log2(l_m / l_box)
    #
    #    We need L' to be a positive integer for all cells.
    #    The subbox side must be a power-of-two fraction of the
    #    full domain for the octree to be consistent. In general
    #    it won't be, so we round l_box up to the next power-of-two
    #    fraction of l_m and pad accordingly.
    # ----------------------------------------------------------
    l_box_raw = 2.0 * hw_m

    # Find the smallest power-of-two fraction of l_m that contains the subbox
    # l_m / 2^k <= l_box_raw  =>  k = floor(log2(l_m / l_box_raw))
    k = int(np.floor(np.log2(l_m / l_box_raw)))
    if k < 0:
        k = 0
    l_box = l_m / (2**k)

    # Level offset: cells at level L in the full domain become level (L - k)
    # in the subbox octree. The minimum level in the subbox should be >= 0.
    level_offset = k
    sub_level_local = (sub_level - level_offset).astype(int)

    if sub_level_local.min() < 0:
        # This means some cells are coarser than the subbox itself.
        # We need a larger subbox. Reduce k.
        min_level_full = int(sub_level.min())
        level_offset = min_level_full
        k = min_level_full
        l_box = l_m / (2**k)

        sub_level_local = (sub_level - level_offset).astype(int)

    # Recompute half-width in case l_box was adjusted
    hw_actual = l_box / 2.0

    # Recenter cell positions relative to the subbox origin
    # Subbox origin = sink position - half the (possibly enlarged) box
    box_origin_x = sink_pos_m[0] - hw_actual
    box_origin_y = sink_pos_m[1] - hw_actual
    box_origin_z = sink_pos_m[2] - hw_actual

    # Cell positions relative to subbox origin
    rel_x = sub_cx - box_origin_x
    rel_y = sub_cy - box_origin_y
    rel_z = sub_cz - box_origin_z

    l_box_cm = l_box * 100.0

    max_level_local = int(sub_level_local.max())
    min_level_local = int(sub_level_local.min())

    print(f"\n--- Subbox for sink {sink_idx} ---")
    print(f"    Sink position (AU): "
          f"({sink_pos_m[0]/au2m:.1f}, {sink_pos_m[1]/au2m:.1f}, {sink_pos_m[2]/au2m:.1f})")
    print(f"    Requested half-width: {box_half_width_au} AU")
    print(f"    Actual subbox side:   {l_box/au2m:.1f} AU "
          f"(rounded to power-of-two fraction of domain)")
    print(f"    Cells in subbox:      {nb_cells_sub}")
    print(f"    Level range (local):  {min_level_local} - {max_level_local}")
    print(f"    Level offset:         {level_offset}")
    print(f"    Dust species:         {nb_species}\n")

    # ----------------------------------------------------------
    # 4) Build octree
    # ----------------------------------------------------------
    tree = OcTree(0.0, 0.0, 0.0, l_box)
    tree.nr_of_cells = nb_cells_sub

    hole_radius2 = (hole_au * au2m)**2

    print("Constructing subbox octree...")
    for i in range(nb_cells_sub):
        # Position relative to subbox origin (which starts at 0,0,0 for the tree)
        c_x = rel_x[i]
        c_y = rel_y[i]
        c_z = rel_z[i]

        cell_dust = sub_dust[i, :].copy()

        # Dig hole around the central sink
        dx = (sub_cx[i] - sink_pos_m[0])
        dy = (sub_cy[i] - sink_pos_m[1])
        dz = (sub_cz[i] - sink_pos_m[2])
        d2 = dx**2 + dy**2 + dz**2
        if d2 <= hole_radius2:
            cell_dust[:] = 0.0

        cell = CellOct(c_x, c_y, c_z, 0, int(sub_level_local[i]))
        cell.data = cell_dust.tolist()

        tree.insertInTree(tree.root, cell, 0)

        if i % 10000 == 0:
            progress = 100.0 * i / nb_cells_sub
            sys.stdout.write(f"    Constructing subbox octree: {progress:.1f}%\r")
            sys.stdout.flush()

    print("    Constructing subbox octree: done\n")

    # Check integrity
    print("    Checking octree integrity...")
    tree.reset_counter()
    check = tree.checkOcTree(tree.root)
    if not check:
        raise RuntimeError(
            f"Octree integrity check failed for sink {sink_idx}!"
        )
    print("    Octree structure: OK\n")

    # ----------------------------------------------------------
    # 5) Write RADMC-3D files
    # ----------------------------------------------------------
    tree.cell_counter = 0
    grid = []
    density = []
    tree._n_species = nb_species
    tree.writeOcTree_radmc(tree.root, grid, density)
    densityarray = np.array(density)

    print(f"    Writing amr_grid.inp ...")
    radmc3d.write.amr_grid(
        output_dir, grid, max_level_local, nb_cells_sub, l_box_cm,
        gridstyle=gridstyle, coordsystem=coordsystem,
    )

    print(f"    Writing dust_density.inp ...")
    radmc3d.write.dust_density(
        output_dir, densityarray, nb_cells_sub, nb_species,
        gridstyle=gridstyle,
    )

    # stars.inp: single central star at the subbox center
    lam = np.logspace(np.log10(lmin), np.log10(lmax), nlam)

    # Star position in cm relative to the subbox center (which is the sink)
    star_pos_cm = np.array([[0.0, 0.0, 0.0]])
    star_mass = np.array([sink_mass_g])
    star_radius = np.array([sink_radius])
    star_teff = np.array([sink_teff])

    print(f"    Writing stars.inp ...")
    radmc3d.write.stars(
        output_dir, 1, star_mass, star_pos_cm, star_radius, lam, star_teff,
    )

    print(f"    Sink {sink_idx} done -> {output_dir}\n")

    return grid, densityarray