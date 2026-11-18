"""
_____________________________________________________________________________________________________________
file name: dustbins
last update: Oct 2026
language: > PYTHON 3.9
short description: Resolve the dust mixtures written in the POLARIS command files against the number of
                   dust density fields in the grid. For multi-fluid RAMSES outputs, a single mixture is
                   used as a composition template and replicated over the dust size bins read from
                   info_mf_XXXXX.txt (one mixture per fluid, in fluid order), with each bin clipped to the
                   grain-size range covered by the tabulated cross sections.
_____________________________________________________________________________________________________________
"""

import copy
from pathlib import Path

import numpy as np

from radprocess import ramses


# ============================================================
#  Size bins
# ============================================================

def bin_edges_from_centres(centres):
    """
    Log-spaced bin edges from the geometric bin centres (RAMSES mf_grain_size_i).

    Interior edges are the geometric means of neighbouring centres; the two
    outer edges are extrapolated with the same log spacing.
    """
    c = np.asarray(centres, dtype=float)
    if c.size < 2:
        raise ValueError("At least two bin centres are needed to derive bin edges.")

    inner = np.sqrt(c[:-1] * c[1:])
    first = c[0] ** 2 / inner[0]
    last = c[-1] ** 2 / inner[-1]
    return np.concatenate([[first], inner, [last]])


def table_grain_sizes(path):
    """
    Grain sizes [m] tabulated in a POLARIS cross-section file (#a_eff block).

    Returns None for refractive-index files (.nk), for which POLARIS computes
    Mie cross sections at any size, or if no #a_eff block is found.
    """
    path = Path(path)
    if path.suffix.lower() == ".nk":
        return None

    with open(path, "r") as f:
        in_block = False
        for line in f:
            stripped = line.strip()
            if stripped.startswith("#a_eff"):
                in_block = True
                continue
            if in_block and stripped and not stripped.startswith("#"):
                return np.array([float(v) for v in stripped.split()])
    return None


def clip_bin(amin, amax, a_eff):
    """
    Clip the size bin [amin, amax] to the tabulated range a_eff.

    A bin lying entirely outside the table is mapped onto the first or last
    tabulated size interval, so it keeps a valid (non-empty) size range and its
    one-to-one correspondence with its dust density field.
    """
    lo, hi = a_eff[0], a_eff[-1]
    if amin >= hi:
        return a_eff[-2], hi
    if amax <= lo:
        return lo, a_eff[1]
    return max(amin, lo), min(amax, hi)


# ============================================================
#  Mixtures
# ============================================================

def expand_dust_mixtures(template, edges):
    """
    Replicate one mixture over the size bins defined by `edges`.

    Parameters
    ----------
    template : dict
        One mixture, {component_id: component_dict}. amin/amax are not used.
    edges : array-like
        Bin edges [m], length n_bins + 1.

    Returns
    -------
    mixtures : dict
        {i_bin: {component_id: component_dict}} with amin/amax set per bin.
    clipped : list of (i_bin, (amin, amax), (amin_clipped, amax_clipped))
        Bins whose size range was clipped to a tabulated cross-section range.
    """
    tables = {}
    mixtures = {}
    clipped = []

    for i in range(len(edges) - 1):
        components = copy.deepcopy(template)
        for component in components.values():
            amin, amax = float(edges[i]), float(edges[i + 1])

            path = component["path"]
            if path not in tables:
                tables[path] = table_grain_sizes(path)
            a_eff = tables[path]

            if a_eff is not None and len(a_eff) > 1:
                new = clip_bin(amin, amax, a_eff)
                clip_reported = clipped and clipped[-1][0] == i
                if not np.allclose(new, (amin, amax), rtol=1e-6) and not clip_reported:
                    clipped.append((i, (amin, amax), new))
                amin, amax = new

            component["amin"] = amin
            component["amax"] = amax
        mixtures[i] = components

    return mixtures, clipped


def _check_sizes_given(dust_mixtures):
    """Raise if a component lacks amin/amax (required when mixtures are used as-is)."""
    for mixture, components in dust_mixtures.items():
        for name, component in components.items():
            missing = [k for k in ("amin", "amax") if k not in component]
            if missing:
                raise ValueError(
                    f"Dust mixture {mixture}, component {name}: missing {', '.join(missing)}. "
                    "Grain sizes are required unless a single mixture is replicated over "
                    "multi-fluid size bins."
                )


def check_mie_nr_sizes(mie_nr_sizes, dust_mixtures, verbose=True):
    """
    Validate the number of Mie grain sizes written as POLARIS <mie_nr_dust_size>.

    Parameters
    ----------
    mie_nr_sizes : int or None
        Number of grain sizes per dust component. None keeps the POLARIS default.
    dust_mixtures : dict
        {mixture_id: {component_id: component_dict}} (see DustConfig).
    verbose : bool
        Print a note if no component is computed with Mie theory.
    """
    if mie_nr_sizes is None:
        return

    if isinstance(mie_nr_sizes, bool) or not isinstance(mie_nr_sizes, (int, np.integer)) or mie_nr_sizes < 2:
        raise ValueError(f"mie_nr_sizes must be an integer >= 2, got {mie_nr_sizes!r}.")

    uses_mie = any(
        Path(component["path"]).suffix.lower() == ".nk"
        for components in dust_mixtures.values()
        for component in components.values()
    )
    if verbose and not uses_mie:
        print(f"NOTE: mie_nr_sizes = {mie_nr_sizes} has no effect: no dust component uses an .nk file; "
              "tabulated cross sections use the grain sizes of their file.")


def resolve_dust_mixtures(dust_mixtures, n_dust, ramses_dir, bin_edges=None, verbose=True):
    """
    Return the dust mixtures to write in the POLARIS command files, one per
    dust density field in the grid.

    * len(dust_mixtures) == n_dust : used as-is (amin/amax required). This is
      the single-species case and any hand-written set of mixtures.
    * len(dust_mixtures) == 1 and n_dust > 1 : the mixture is a composition
      template replicated over the n_dust size bins, taken from `bin_edges`
      or, if None, from info_mf_XXXXX.txt in `ramses_dir`.

    Parameters
    ----------
    dust_mixtures : dict
        {mixture_id: {component_id: component_dict}} (see DustConfig).
    n_dust : int
        Number of dust density fields in the grid.
    ramses_dir : str or Path
        RAMSES output directory (for info_mf_XXXXX.txt and namelist.txt).
    bin_edges : list or None
        Explicit bin edges [m], length n_dust + 1.
    verbose : bool
        Print the resulting size bins and any clipping.

    Returns
    -------
    mixtures : dict
    """
    n_mix = len(dust_mixtures)

    if n_mix == n_dust:
        _check_sizes_given(dust_mixtures)
        return dust_mixtures

    if n_mix != 1:
        raise ValueError(
            f"{n_mix} dust mixtures are defined but the grid holds {n_dust} dust density "
            f"fields. Define either a single mixture (replicated over the size bins) or "
            f"exactly {n_dust} mixtures."
        )

    # ---- one template, n_dust bins ----
    grain_density_kgm3 = None
    if bin_edges is not None:
        edges = np.asarray(bin_edges, dtype=float)
        source = "configparams.dust.bin_edges"
    else:
        mf = ramses.read.mf_grain_sizes(ramses_dir)
        if mf is None:
            raise ValueError(
                f"The grid holds {n_dust} dust density fields but only one dust mixture is "
                f"defined, and no info_mf_*.txt was found in {ramses_dir} to derive the size "
                f"bins. Set configparams.dust.bin_edges ({n_dust + 1} edges in m) or define "
                f"{n_dust} mixtures."
            )
        if len(mf["sizes_cm"]) != n_dust:
            raise ValueError(
                f"{Path(mf['filename']).name} lists {len(mf['sizes_cm'])} grain sizes but the "
                f"grid holds {n_dust} dust density fields."
            )
        edges = bin_edges_from_centres(mf["sizes_cm"] * 1e-2)  # cm -> m
        source = Path(mf["filename"]).name
        grain_density_kgm3 = np.unique(mf["densities_gcc"][np.isfinite(mf["densities_gcc"])]) * 1e3

        # Cross-check with the namelist size range, when available
        size_min = ramses.read.namelist_value(ramses_dir, "mf_size_min")
        size_max = ramses.read.namelist_value(ramses_dir, "mf_size_max")
        if size_min is not None and size_max is not None:
            if not np.allclose(edges[[0, -1]], np.array([size_min, size_max]) * 1e-2, rtol=1e-3):
                print(f"WARNING: size bins from {source} span {edges[0]:.3e}-{edges[-1]:.3e} m, "
                      f"but namelist.txt has mf_size_min/max = {size_min * 1e-2:.3e}/"
                      f"{size_max * 1e-2:.3e} m.")

    if len(edges) != n_dust + 1:
        raise ValueError(f"{len(edges)} bin edges given for {n_dust} dust density fields "
                         f"(expected {n_dust + 1}).")
    if np.any(np.diff(edges) <= 0):
        raise ValueError("Dust bin edges must be strictly increasing.")

    template = next(iter(dust_mixtures.values()))
    mixtures, clipped = expand_dust_mixtures(template, edges)

    if verbose:
        print(f"\nDust size bins: 1 mixture replicated over {n_dust} bins ({source})")
        print(f"  range: {edges[0]:.3e} - {edges[-1]:.3e} m")

        if any("amin" in c or "amax" in c for c in template.values()):
            print("  NOTE: amin/amax of the template mixture are ignored; the bin edges are used.")

        if grain_density_kgm3 is not None:
            for name, component in template.items():
                if component.get("density", 0) == 0:
                    print(f"  NOTE: component {name} uses the material density of its dust file; "
                          f"RAMSES used {', '.join(f'{d:g}' for d in grain_density_kgm3)} kg/m^3.")
                elif not np.allclose(component["density"], grain_density_kgm3, rtol=1e-3):
                    print(f"  WARNING: component {name} density = {component['density']:g} kg/m^3, "
                          f"RAMSES used {', '.join(f'{d:g}' for d in grain_density_kgm3)} kg/m^3.")

        if clipped:
            print(f"  WARNING: {len(clipped)} bin(s) clipped to the tabulated grain-size range:")
            for i, (amin, amax), (cmin, cmax) in clipped:
                print(f"    bin {i + 1:3d}: [{amin:.3e}, {amax:.3e}] -> [{cmin:.3e}, {cmax:.3e}] m")
        print()

    return mixtures
