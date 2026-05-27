#!/usr/bin/env python
# -*- coding: utf-8 -*-
import math

import numpy as np
import matplotlib.pyplot as plt

from ..constants.constants import au2cm, M_sun, R_sun


def sink_columns(sink_data, x_col, y_col):
    """
    sink_data: list of dicts from sink_info
    x_col, y_col: strings
    """
    x = [float(row[x_col]) for row in sink_data]
    y = [float(row[y_col]) for row in sink_data]

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(x, y)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title(f"{y_col} vs {x_col}")
    return fig



def ramses_histogram(root, y_field):
    """
    Scatter-density histogram of dust density vs chosen RAMSES field.

    Parameters
    ----------
    root : zarr.Group
        RAMSES AMR Zarr root
    y_field : str
        Field name for Y axis (e.g. gas_massdensity, Tgas, etc.)

    Returns
    -------
    matplotlib.figure.Figure
    """

    # ----------------------------
    # Load Y
    # ----------------------------
    if y_field not in root:
        raise RuntimeError(f"Field '{y_field}' not found in Zarr.")

    Y = root[y_field][:]

    # ----------------------------
    # Load dust
    # ----------------------------
    if "dust_massdensities" in root:
        # multi species: (cells, nspecies)
        dust = root["dust_massdensities"][:]
        nspecies = dust.shape[1]
    elif "dust_massdensity" in root:
        # single species: (cells,)
        dust = root["dust_massdensity"][:][:, None]
        nspecies = 1
    else:
        raise RuntimeError("No dust_massdensity(s) found in Zarr.")

    # ----------------------------
    # Log safety
    # ----------------------------
    dust = np.clip(dust, 1e-99, None)
    Y = np.clip(Y, 1e-99, None)

    logdust = np.log10(dust)
    logY = np.log10(Y)

    # ----------------------------
    # Figure layout
    # ----------------------------
    ncols = 4
    nrows = math.ceil(nspecies / ncols)

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4 * ncols, 3 * nrows),
        squeeze=False,
    )

    # ----------------------------
    # Plot each species
    # ----------------------------
    for i in range(nspecies):
        r = i // ncols
        c = i % ncols
        ax = axes[r, c]

        x = logdust[:, i]
        y = logY

        hb = ax.hexbin(
            x,
            y,
            gridsize=150,
            bins="log",
            mincnt=1,
        )

        ax.set_xlabel(r"log$_{10}(\rho_\mathrm{dust})$")
        ax.set_ylabel(f"log$_{{10}}$({y_field})")
        ax.set_title(f"Species {i+1}")

    # turn off unused panels
    for j in range(nspecies, nrows * ncols):
        r = j // ncols
        c = j % ncols
        axes[r, c].axis("off")

    plt.tight_layout()

    return fig


def density_midplane(root, plane="xy", field="gas", npix=512, log=True,
                     cmap="inferno", unit_label=None, vmin=None, vmax=None):
    """
    Plot a 2D midplane density map from the Zarr grid data.

    This works for checking both the RADMC-3D and POLARIS grids,
    since both are built from the same Zarr data.

    The function selects cells closest to the midplane (within one
    cell width of the plane center), projects them onto a 2D pixel
    grid, and displays the result.

    Parameters
    ----------
    root : zarr.Group
        RAMSES AMR Zarr root (from pipe.get_amr_root()).
    plane : str
        Projection plane: "xy", "xz", or "yz".
    field : str
        Which density to plot:
        - "gas": gas mass density
        - "dust": total dust density (sum of all species)
        - "dust_0", "dust_1", ...: individual dust species
        - "all_dust": one subplot per dust species
    npix : int
        Number of pixels per side for the 2D map.
    log : bool
        If True, plot log10 of the density.
    cmap : str
        Matplotlib colormap name.
    unit_label : str or None
        Label for the colorbar. If None, auto-generated.
    vmin, vmax : float or None
        Colorbar limits (applied to log10 values if log=True).

    Returns
    -------
    matplotlib.figure.Figure
    """
    from ..constants.constants import au2m

    # Load coordinates and cell sizes
    x = np.array(root["x"])
    y = np.array(root["y"])
    z = np.array(root["z"])
    dx = np.array(root["dx"])
    l_m = root.attrs.get("l_m")

    # Center coordinates
    cx = x - 0.5 * l_m
    cy = y - 0.5 * l_m
    cz = z - 0.5 * l_m

    # Convert to AU for display
    cx_au = cx / au2m
    cy_au = cy / au2m
    cz_au = cz / au2m
    dx_au = dx / au2m
    half_box_au = 0.5 * l_m / au2m

    # Load densities
    if "dust_massdensities" in root:
        dust_all = np.array(root["dust_massdensities"])
        n_dust = dust_all.shape[1]
    elif "dust_massdensity" in root:
        dust_all = np.array(root["dust_massdensity"])
        if dust_all.ndim == 1:
            dust_all = dust_all[:, np.newaxis]
        n_dust = 1
    else:
        dust_all = None
        n_dust = 0

    gas = np.array(root["gas_massdensity"]) if "gas_massdensity" in root else None

    # Determine which fields to plot
    if field == "all_dust":
        if dust_all is None:
            raise RuntimeError("No dust density found in Zarr.")
        fields_to_plot = [(f"Dust species {i+1}", dust_all[:, i]) for i in range(n_dust)]
    elif field == "gas":
        if gas is None:
            raise RuntimeError("No gas_massdensity found in Zarr.")
        fields_to_plot = [("Gas density", gas)]
    elif field == "dust":
        if dust_all is None:
            raise RuntimeError("No dust density found in Zarr.")
        fields_to_plot = [("Total dust density", dust_all.sum(axis=1))]
    elif field.startswith("dust_"):
        idx = int(field.split("_")[1])
        if dust_all is None or idx >= n_dust:
            raise RuntimeError(f"Dust species {idx} not found (have {n_dust}).")
        fields_to_plot = [(f"Dust species {idx+1}", dust_all[:, idx])]
    else:
        # Try to load it as a generic field
        if field in root:
            fields_to_plot = [(field, np.array(root[field]))]
        else:
            raise RuntimeError(f"Field '{field}' not found in Zarr.")

    # Select midplane axis
    if plane == "xy":
        coord_perp = cz
        c1, c2 = cx_au, cy_au
        label1, label2 = "x [AU]", "y [AU]"
    elif plane == "xz":
        coord_perp = cy
        c1, c2 = cx_au, cz_au
        label1, label2 = "x [AU]", "z [AU]"
    elif plane == "yz":
        coord_perp = cx
        c1, c2 = cy_au, cz_au
        label1, label2 = "y [AU]", "z [AU]"
    else:
        raise ValueError(f"Unknown plane '{plane}'. Use 'xy', 'xz', or 'yz'.")

    # Select cells near the midplane (within one cell width of center)
    midplane_mask = np.abs(coord_perp) <= dx

    if midplane_mask.sum() == 0:
        raise RuntimeError(f"No cells found near the {plane} midplane.")

    # Figure layout
    n_plots = len(fields_to_plot)
    ncols = min(n_plots, 4)
    nrows = math.ceil(n_plots / ncols)

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(5 * ncols, 4.5 * nrows),
        squeeze=False,
    )

    for idx, (title, density_field) in enumerate(fields_to_plot):
        r = idx // ncols
        c = idx % ncols
        ax = axes[r, c]

        # Build 2D histogram (density-weighted)
        sel_c1 = c1[midplane_mask]
        sel_c2 = c2[midplane_mask]
        sel_dens = density_field[midplane_mask]
        sel_dx = dx_au[midplane_mask]

        # Weight by cell area (dx^2) to get proper averaging
        weights = sel_dens * sel_dx**2

        extent = [-half_box_au, half_box_au, -half_box_au, half_box_au]

        img, xedges, yedges = np.histogram2d(
            sel_c1, sel_c2,
            bins=npix,
            range=[[-half_box_au, half_box_au], [-half_box_au, half_box_au]],
            weights=weights,
        )

        # Normalize by pixel area contribution
        area, _, _ = np.histogram2d(
            sel_c1, sel_c2,
            bins=npix,
            range=[[-half_box_au, half_box_au], [-half_box_au, half_box_au]],
            weights=sel_dx**2,
        )

        with np.errstate(invalid="ignore", divide="ignore"):
            img = np.where(area > 0, img / area, np.nan)

        if log:
            with np.errstate(invalid="ignore", divide="ignore"):
                img = np.log10(img)
            if unit_label is None:
                clabel = r"log$_{10}$($\rho$) [g/cm$^3$]"
            else:
                clabel = unit_label
        else:
            if unit_label is None:
                clabel = r"$\rho$ [g/cm$^3$]"
            else:
                clabel = unit_label

        im = ax.imshow(
            img.T,
            origin="lower",
            extent=extent,
            aspect="equal",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )

        ax.set_xlabel(label1)
        ax.set_ylabel(label2)
        ax.set_title(title)
        fig.colorbar(im, ax=ax, label=clabel, shrink=0.85)

    # Turn off unused panels
    for j in range(n_plots, nrows * ncols):
        r = j // ncols
        c = j % ncols
        axes[r, c].axis("off")

    plt.tight_layout()
    return fig


# ============================================================
#  Subbox dust density mosaic
# ============================================================

def _read_dust_density_subbox(filepath):
    """
    Read a RADMC-3D dust_density_subbox.out file.

    Returns
    -------
    density : np.ndarray of shape (nx, ny, nz)
        Dust density in g/cm^3 (Fortran ordering).
    bounds : tuple
        (xmin, xmax, ymin, ymax, zmin, zmax) in cm.
    """
    with open(filepath, "r") as f:
        lines = f.readlines()

    lines = [l.strip() for l in lines if l.strip() != ""]

    # line 0: format number
    # line 1: nx ny nz nrspec
    nx, ny, nz, ndust = map(int, lines[1].split())

    # line 2: xmin xmax ymin ymax zmin zmax
    xmin, xmax, ymin, ymax, zmin, zmax = map(float, lines[2].split())
    bounds = (xmin, xmax, ymin, ymax, zmin, zmax)

    # Data starts at line 5
    data = np.array([float(x) for x in lines[5:]])

    expected_size = nx * ny * nz
    if data.size != expected_size:
        raise ValueError(
            f"Expected {expected_size} values but found {data.size}"
        )

    density = data.reshape((nx, ny, nz), order="F")
    return density, bounds


def _project_density(density, axis="z"):
    """
    Column density projection along chosen axis.

    Parameters
    ----------
    density : np.ndarray
        3D density cube with shape (nx, ny, nz).
    axis : str
        'x', 'y', or 'z'.

    Returns
    -------
    map2d : np.ndarray
    """
    axis_dict = {"x": 0, "y": 1, "z": 2}
    if axis not in axis_dict:
        raise ValueError("axis must be 'x', 'y', or 'z'")
    return np.sum(density, axis=axis_dict[axis])


def subbox_density_mosaic(
    subbox_base_dir,
    axis="z",
    n_sinks=10,
    vmin=None,
    vmax=None,
    log=True,
    cmap="inferno",
    figsize=None,
    fontsize=16,
    sink_folders=None,
):
    """
    Plot a mosaic of projected dust column density maps from RADMC-3D
    subbox regrid output, with a single shared colorbar.

    Parameters
    ----------
    subbox_base_dir : str or Path
        Path to radmc3d/subboxes/ containing sink_0XXX/ folders.
    axis : str
        Projection axis: "x", "y", or "z" (default "z").
    n_sinks : int
        Maximum number of sinks to plot (default 10).
    vmin, vmax : float or None
        Colorbar limits (applied to log10 values if log=True).
    log : bool
        If True, plot log10 of the column density (default True).
    cmap : str
        Matplotlib colormap (default "inferno").
    figsize : tuple or None
        Figure size. If None, auto-computed.
    fontsize : int
        Font size for titles and colorbar label (default 16).
    sink_folders : list of str or None
        If provided, only plot these specific folders.
        Overrides n_sinks.

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    from pathlib import Path

    subbox_base_dir = Path(subbox_base_dir)

    # Find sink folders with dust_density_subbox.out
    if sink_folders is not None:
        folders = []
        for name in sink_folders:
            d = subbox_base_dir / name
            if (d / "dust_density_subbox.out").exists():
                folders.append(d)
            else:
                print(f"WARNING: {d / 'dust_density_subbox.out'} not found, skipping.")
    else:
        folders = sorted([
            d for d in subbox_base_dir.iterdir()
            if d.is_dir() and d.name.startswith("sink_")
            and (d / "dust_density_subbox.out").exists()
        ])[:n_sinks]

    if not folders:
        raise FileNotFoundError(
            f"No dust_density_subbox.out files found in {subbox_base_dir}."
        )

    n_plots = len(folders)
    ncols = min(n_plots, 5)
    nrows = math.ceil(n_plots / ncols)

    if figsize is None:
        figsize = (3.5 * ncols, 3.5 * nrows)

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=figsize,
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes).ravel()

    images = []

    for i, folder in enumerate(folders):
        filepath = folder / "dust_density_subbox.out"
        density, bounds = _read_dust_density_subbox(filepath)

        img = _project_density(density, axis=axis)

        if log:
            with np.errstate(invalid="ignore", divide="ignore"):
                img = np.log10(img + 1e-99)

        # Determine the 2D extent from the 3D bounds
        # bounds = (xmin, xmax, ymin, ymax, zmin, zmax) in cm
        xmin, xmax, ymin, ymax, zmin, zmax = bounds
        au2cm_local = 1.495978707e13
        if axis == "z":
            ext = [xmin/au2cm_local, xmax/au2cm_local,
                   ymin/au2cm_local, ymax/au2cm_local]
        elif axis == "y":
            ext = [xmin/au2cm_local, xmax/au2cm_local,
                   zmin/au2cm_local, zmax/au2cm_local]
        elif axis == "x":
            ext = [ymin/au2cm_local, ymax/au2cm_local,
                   zmin/au2cm_local, zmax/au2cm_local]
        else:
            ext = None

        im = axes[i].imshow(
            img.T,
            origin="lower",
            extent=ext,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
        images.append(im)

        # Mark the sink position at the grid center (0,0,0)
        axes[i].plot(0, 0, "w+", ms=8, mew=1.5)

        axes[i].set_title(folder.name, fontsize=fontsize, fontweight="bold")
        axes[i].tick_params(labelsize=fontsize - 4)

    # Hide unused axes
    for j in range(n_plots, len(axes)):
        axes[j].axis("off")

    # Single shared colorbar
    if log:
        clabel = r"$\log_{10}(\Sigma_{\rm dust})$  [g/cm$^2$]"
    else:
        clabel = r"$\Sigma_{\rm dust}$  [g/cm$^2$]"

    cbar = fig.colorbar(
        images[0],
        ax=axes.tolist(),
        shrink=0.85,
        pad=0.02,
    )
    cbar.set_label(clabel, fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize - 2)

    return fig