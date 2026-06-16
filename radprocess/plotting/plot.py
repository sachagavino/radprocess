#!/usr/bin/env python
# -*- coding: utf-8 -*-
import math
from pathlib import Path

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

        # Determine the 2D extent from the 3D bounds, shifted to
        # sink-centered coordinates so the sink is at (0,0) in all subplots.
        # bounds = (xmin, xmax, ymin, ymax, zmin, zmax) in cm
        xmin, xmax, ymin, ymax, zmin, zmax = bounds
        au2cm_local = 1.495978707e13

        # Read the sink offset to shift coordinates
        offset_file = folder / "sink_offset.txt"
        if offset_file.exists():
            offset = np.loadtxt(offset_file)
            ox, oy, oz = offset  # in cm (RADMC-3D subboxes)
        else:
            ox, oy, oz = 0.0, 0.0, 0.0

        if axis == "z":
            ext = [(xmin - ox)/au2cm_local, (xmax - ox)/au2cm_local,
                   (ymin - oy)/au2cm_local, (ymax - oy)/au2cm_local]
        elif axis == "y":
            ext = [(xmin - ox)/au2cm_local, (xmax - ox)/au2cm_local,
                   (zmin - oz)/au2cm_local, (zmax - oz)/au2cm_local]
        elif axis == "x":
            ext = [(ymin - oy)/au2cm_local, (ymax - oy)/au2cm_local,
                   (zmin - oz)/au2cm_local, (zmax - oz)/au2cm_local]
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

        # Sink is at (0,0) in sink-centered coordinates
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


# ============================================================
#  Unified subbox mosaic (density, intensity, temperature)
# ============================================================

def _read_polaris_fits(filepath, stokes=0):
    """
    Read a POLARIS detector FITS file and return Stokes I (or Q, U, V).

    Parameters
    ----------
    filepath : str or Path
        Path to polaris_detector_nr*.fits.gz
    stokes : int
        Stokes index: 0=I, 1=Q, 2=U, 3=V, 4=tau, 5=column_density

    Returns
    -------
    image : np.ndarray (ny, nx)
    extent_m : list [xmin, xmax, ymin, ymax] in metres
    """
    from astropy.io import fits as pyfits

    hdul = pyfits.open(filepath)
    data = hdul[0].data           # shape (6, 1, ny, nx)
    header = hdul[0].header
    hdul.close()

    image = data[stokes, 0, :, :]  # (ny, nx)

    cdelt1 = header["CDELT1"]      # pixel size in m
    cdelt2 = header["CDELT2"]
    nx = header["NAXIS1"]
    ny = header["NAXIS2"]

    half_x = nx * cdelt1 / 2.0
    half_y = ny * cdelt2 / 2.0
    extent_m = [-half_x, half_x, -half_y, half_y]

    return image, extent_m


def _find_sink_offset(sink_name, pipeline_output):
    """
    Find and read sink_offset.txt for a given sink.
    Returns (ox, oy, oz) in metres.
    """
    pipeline_output = Path(pipeline_output)
    au2m_local = 1.4959787070e11

    # Try POLARIS subbox (offset in metres)
    polaris_offset = pipeline_output / "polaris" / "subboxes" / sink_name / "sink_offset.txt"
    if polaris_offset.exists():
        return np.loadtxt(polaris_offset)

    # Try RADMC-3D subbox (offset in cm → convert to m)
    radmc_offset = pipeline_output / "radmc3d" / "subboxes" / sink_name / "sink_offset.txt"
    if radmc_offset.exists():
        return np.loadtxt(radmc_offset) / 100.0  # cm → m

    return np.array([0.0, 0.0, 0.0])


def subbox_mosaic(
    pipeline_output,
    quantity="intensity",
    view="xy",
    wavelength_idx=0,
    stokes=0,
    n_sinks=10,
    vmin=None,
    vmax=None,
    log=True,
    cmap="inferno",
    figsize=None,
    fontsize=16,
    sink_folders=None,
    label="whole",
):
    """
    Plot a mosaic of subbox images: intensity (POLARIS), density (RADMC-3D
    subbox_regrid), or optical depth, centered on each sink.

    Parameters
    ----------
    pipeline_output : str or Path
        Pipeline output directory (e.g., test1/).
    quantity : str
        What to plot:
        - "intensity": Stokes I from POLARIS FITS (Jy/px)
        - "density": column density from dust_density_subbox.out (g/cm^2)
        - "tau": optical depth from POLARIS FITS
        - "column_density": column density from POLARIS FITS (m^-2)
    view : str
        Viewing direction: "xy", "xz", or "yz".
    wavelength_idx : int
        Which wavelength to plot (0 = first). Only for POLARIS quantities.
    stokes : int
        Stokes parameter: 0=I, 1=Q, 2=U, 3=V. Only for "intensity".
    n_sinks : int
        Maximum number of sinks to plot.
    vmin, vmax : float or None
        Colorbar limits (applied to log10 values if log=True).
    log : bool
        If True, plot log10 of the quantity.
    cmap : str
        Matplotlib colormap.
    figsize : tuple or None
        Figure size. If None, auto-computed.
    fontsize : int
        Font size for titles and labels.
    sink_folders : list of str or None
        If provided, only plot these specific sink folders.
    label : str
        Subfolder label for POLARIS images (default "whole").

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    pipeline_output = Path(pipeline_output)
    au2m_local = 1.4959787070e11
    au2cm_local = 1.495978707e13

    # Map quantity to POLARIS FITS index
    fits_index_map = {"intensity": 0, "tau": 4, "column_density": 5}

    if quantity == "density":
        # Use RADMC-3D subbox_regrid output
        subbox_base = pipeline_output / "radmc3d" / "subboxes"
        file_check = "dust_density_subbox.out"
    elif quantity in fits_index_map:
        # Use POLARIS FITS output
        subbox_base = pipeline_output / "images" / "subboxes"
        file_check = None  # check for FITS files below
        if quantity == "intensity":
            fits_stokes = stokes
        else:
            fits_stokes = fits_index_map[quantity]
    else:
        raise ValueError(f"quantity must be 'intensity', 'density', 'tau', "
                         f"or 'column_density', got '{quantity}'")

    # Find sink folders
    if sink_folders is not None:
        folder_names = sink_folders
    else:
        folder_names = sorted([
            d.name for d in subbox_base.iterdir()
            if d.is_dir() and d.name.startswith("sink_")
        ])[:n_sinks]

    # Validate folders and find data
    folders = []
    for name in folder_names:
        if quantity == "density":
            data_path = subbox_base / name / file_check
            if data_path.exists():
                folders.append(name)
            else:
                print(f"WARNING: {data_path} not found, skipping {name}")
        else:
            fits_dir = subbox_base / name / label / view / "data"
            fits_pattern = f"polaris_detector_nr{wavelength_idx+1:04d}.fits*"
            matches = sorted(fits_dir.glob(fits_pattern)) if fits_dir.exists() else []
            if matches:
                folders.append(name)
            else:
                print(f"WARNING: no FITS for {name}/{label}/{view}/, skipping")

    if not folders:
        raise FileNotFoundError(
            f"No {quantity} data found in {subbox_base}."
        )

    n_plots = len(folders)
    ncols = min(n_plots, 5)
    nrows = math.ceil(n_plots / ncols)

    if figsize is None:
        figsize = (3.5 * ncols, 3.5 * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize,
                              constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()

    images = []

    for i, sink_name in enumerate(folders):
        # Read sink offset (in metres)
        offset_m = _find_sink_offset(sink_name, pipeline_output)
        ox_m, oy_m, oz_m = offset_m

        if quantity == "density":
            # --- RADMC-3D density subbox ---
            filepath = subbox_base / sink_name / "dust_density_subbox.out"
            density, bounds = _read_dust_density_subbox(filepath)
            axis_for_proj = {"xy": "z", "xz": "y", "yz": "x"}[view]
            img = _project_density(density, axis=axis_for_proj)

            xmin, xmax, ymin, ymax, zmin, zmax = bounds
            # Offset in cm for RADMC-3D
            ox_cm, oy_cm, oz_cm = offset_m * 100.0
            if view == "xy":
                ext = [(xmin - ox_cm)/au2cm_local, (xmax - ox_cm)/au2cm_local,
                       (ymin - oy_cm)/au2cm_local, (ymax - oy_cm)/au2cm_local]
            elif view == "xz":
                ext = [(xmin - ox_cm)/au2cm_local, (xmax - ox_cm)/au2cm_local,
                       (zmin - oz_cm)/au2cm_local, (zmax - oz_cm)/au2cm_local]
            elif view == "yz":
                ext = [(ymin - oy_cm)/au2cm_local, (ymax - oy_cm)/au2cm_local,
                       (zmin - oz_cm)/au2cm_local, (zmax - oz_cm)/au2cm_local]

        else:
            # --- POLARIS FITS ---
            fits_dir = subbox_base / sink_name / label / view / "data"
            fits_file = sorted(fits_dir.glob(
                f"polaris_detector_nr{wavelength_idx+1:04d}.fits*"
            ))[0]

            img, extent_m = _read_polaris_fits(fits_file, stokes=fits_stokes)

            # The POLARIS image is already centered on the sink
            # (via detector_shift_m in the <detector_dust> line).
            # The extent from the FITS header is symmetric around
            # the image center = the sink. No offset subtraction needed.
            ext = [extent_m[0]/au2m_local, extent_m[1]/au2m_local,
                   extent_m[2]/au2m_local, extent_m[3]/au2m_local]

        if log:
            with np.errstate(invalid="ignore", divide="ignore"):
                img = np.log10(img + 1e-99)

        im = axes[i].imshow(
            img,
            origin="lower",
            extent=ext,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
        images.append(im)

        axes[i].plot(0, 0, "w+", ms=8, mew=1.5)
        axes[i].set_title(sink_name, fontsize=fontsize, fontweight="bold")
        axes[i].tick_params(labelsize=fontsize - 4)

    for j in range(n_plots, len(axes)):
        axes[j].axis("off")

    # Colorbar
    stokes_labels = {0: "I", 1: "Q", 2: "U", 3: "V"}
    clabels = {
        "intensity": f"Stokes {stokes_labels.get(stokes, stokes)} [Jy/px]",
        "density": r"$\Sigma_{\rm dust}$  [g/cm$^2$]",
        "tau": r"$\tau$",
        "column_density": r"$N$  [m$^{-2}$]",
    }
    clabel = clabels[quantity]
    if log:
        clabel = r"$\log_{10}$(" + clabel + ")"

    cbar = fig.colorbar(images[0], ax=axes.tolist(), shrink=0.85, pad=0.02)
    cbar.set_label(clabel, fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize - 2)

    return fig


# ============================================================
#  Bolometric luminosity and temperature from POLARIS images
# ============================================================

def compute_tbol_lbol(
    image_dir,
    view="xy",
    distance_pc=140.0,
    aperture_au=None,
    label="whole",
):
    """
    Compute T_bol and L_bol from multi-wavelength POLARIS images.

    T_bol = 1.25e-11 × <ν>  [K]
    where <ν> = ∫ ν F_ν dν / ∫ F_ν dν   (Myers & Ladd 1993)

    L_bol = 4π d² ∫ F_ν dν

    Parameters
    ----------
    image_dir : str or Path
        Path to the sink's image directory, e.g.,
        test1/images/subboxes/sink_0026/
    view : str
        Viewing direction: "xy", "xz", or "yz".
    distance_pc : float
        Source distance in parsecs.
    aperture_au : float or None
        Aperture radius in AU. If None, integrate over the full image.
    label : str
        Subdirectory label (default "whole").

    Returns
    -------
    result : dict
        "T_bol": bolometric temperature [K]
        "L_bol": bolometric luminosity [L_sun]
        "wavelengths_um": wavelength array [µm]
        "flux_Jy": integrated flux at each wavelength [Jy]
        "aperture_au": aperture used [AU]
    """
    from astropy.io import fits as pyfits

    image_dir = Path(image_dir)
    data_dir = image_dir / label / view / "data"

    # Find all detector FITS files (one per wavelength), skip SED files
    fits_files = sorted(data_dir.glob("polaris_detector_nr*.fits*"))
    fits_files = [f for f in fits_files if "_sed" not in f.name]

    if not fits_files:
        raise FileNotFoundError(
            f"No polaris_detector_nr*.fits files in {data_dir}"
        )

    wavelengths_m = []
    fluxes_Jy = []

    for fits_file in fits_files:
        hdul = pyfits.open(fits_file)
        header = hdul[0].header
        data = hdul[0].data  # (6, 1, ny, nx)
        hdul.close()

        stokes_I = data[0, 0, :, :]  # Jy/px
        ny, nx = stokes_I.shape

        # Read wavelength from header
        wave_m = None
        for key in ["WAVELENGTH1", "WAVELENGTH", "HIERARCH WAVELENGTH1"]:
            if key in header:
                wave_m = float(header[key])
                break
        if wave_m is None:
            raise ValueError(f"No wavelength found in {fits_file.name}")

        wavelengths_m.append(wave_m)

        # Pixel size for aperture masking
        cdelt_m = header["CDELT1"]  # pixel size in metres
        au2m_local = 1.4959787070e11

        if aperture_au is not None:
            cy_pix, cx_pix = ny / 2.0, nx / 2.0
            yy, xx = np.mgrid[0:ny, 0:nx]
            r_au = np.sqrt(
                ((xx - cx_pix) * cdelt_m)**2 +
                ((yy - cy_pix) * cdelt_m)**2
            ) / au2m_local
            mask = r_au <= aperture_au
            total_flux_Jy = np.sum(stokes_I[mask])
        else:
            total_flux_Jy = np.sum(stokes_I)

        fluxes_Jy.append(total_flux_Jy)

    wavelengths_m = np.array(wavelengths_m)
    fluxes_Jy = np.array(fluxes_Jy)

    # Sort by wavelength
    order = np.argsort(wavelengths_m)
    wavelengths_m = wavelengths_m[order]
    fluxes_Jy = fluxes_Jy[order]

    # Convert to frequency (Hz)
    c_light = 2.99792458e8  # m/s
    freqs_Hz = c_light / wavelengths_m

    # Convert Jy to SI: 1 Jy = 1e-26 W/m²/Hz
    fluxes_SI = fluxes_Jy * 1e-26

    # Integrate in frequency space (flip to increasing ν)
    freqs_sorted = freqs_Hz[::-1]
    fluxes_sorted = fluxes_SI[::-1]

    # <ν> = ∫ ν F_ν dν / ∫ F_ν dν
    integral_nu_Fnu = np.trapz(freqs_sorted * fluxes_sorted, freqs_sorted)
    integral_Fnu = np.trapz(fluxes_sorted, freqs_sorted)

    if integral_Fnu > 0:
        mean_freq = integral_nu_Fnu / integral_Fnu
        T_bol = 1.25e-11 * mean_freq
    else:
        mean_freq = 0.0
        T_bol = 0.0

    # L_bol = 4π d² ∫ F_ν dν
    pc2m_local = 3.085677581e16
    L_sun_local = 3.828e26  # W
    distance_m = distance_pc * pc2m_local
    L_bol_W = 4.0 * np.pi * distance_m**2 * integral_Fnu
    L_bol_Lsun = L_bol_W / L_sun_local

    wavelengths_um = wavelengths_m * 1e6

    result = {
        "T_bol": T_bol,
        "L_bol": L_bol_Lsun,
        "wavelengths_um": wavelengths_um,
        "flux_Jy": fluxes_Jy,
        "aperture_au": aperture_au if aperture_au is not None else "full",
        "mean_freq_Hz": mean_freq,
        "integral_Fnu_SI": integral_Fnu,
    }

    print(f"  Wavelength range: {wavelengths_um.min():.1f} - "
          f"{wavelengths_um.max():.0f} um ({len(wavelengths_um)} points)")
    print(f"  Aperture: {result['aperture_au']} AU")
    print(f"  L_bol  = {L_bol_Lsun:.3f} L_sun")
    print(f"  T_bol  = {T_bol:.1f} K")

    return result


def compute_all_tbol_lbol(
    pipeline_output,
    ramses_path,
    view="xy",
    distance_pc=140.0,
    aperture_au=None,
    label="whole",
    sink_folders=None,
):
    """
    Compute T_bol and L_bol for all sinks, and compare with RAMSES values.

    Parameters
    ----------
    pipeline_output : str or Path
        Pipeline output directory (e.g., test1/).
    ramses_path : str or Path
        RAMSES output directory containing the sink CSV file.
    view : str
        Viewing direction: "xy", "xz", or "yz".
    distance_pc : float
        Source distance in parsecs.
    aperture_au : float or None
        Aperture radius in AU. If None, integrate over the full image.
    label : str
        Subdirectory label (default "whole").
    sink_folders : list of str or None
        If provided, only process these sinks.

    Returns
    -------
    table : dict
        Arrays keyed by: "sink_id", "T_bol", "L_bol_obs", "L_bol_true",
        "L_acc", "L_int", "L_bol_ratio" (= L_bol_obs / L_bol_true).
    """
    from radprocess.ramses import read as ramses_read

    pipeline_output = Path(pipeline_output)
    image_base = pipeline_output / "images" / "subboxes"

    # Find sink folders
    if sink_folders is not None:
        folder_names = sink_folders
    else:
        folder_names = sorted([
            d.name for d in image_base.iterdir()
            if d.is_dir() and d.name.startswith("sink_")
        ])

    # Read RAMSES sink data
    sinks = ramses_read.sink_info(str(ramses_path))
    cols = sinks.columns
    id_col = cols.index("Id")
    acclum_col = cols.index("acc_lum[Lsol]")
    intlum_col = cols.index("int_lum[Lsol]")

    # Build lookup: sink_id -> row index
    sink_lookup = {}
    for row_idx in range(sinks.num_sinks):
        sid = int(sinks.data[row_idx, id_col])
        sink_lookup[sid] = row_idx

    # Process each sink
    sink_ids = []
    T_bols = []
    L_bol_obs_list = []
    L_acc_list = []
    L_int_list = []
    L_bol_true_list = []

    for name in folder_names:
        sink_id = int(name.replace("sink_", ""))
        image_dir = image_base / name

        # Check if FITS files exist for this view
        data_dir = image_dir / label / view / "data"
        if not data_dir.exists():
            print(f"  {name}: no {view} data, skipping")
            continue
        fits_files = [f for f in data_dir.glob("polaris_detector_nr*.fits*")
                      if "_sed" not in f.name]
        if not fits_files:
            print(f"  {name}: no FITS files, skipping")
            continue

        # Check if this sink exists in the RAMSES data
        if sink_id not in sink_lookup:
            print(f"  {name}: sink ID {sink_id} not in RAMSES data, skipping")
            continue

        print(f"\n  === {name} (ID={sink_id}) ===")

        try:
            bol = compute_tbol_lbol(
                image_dir, view=view, distance_pc=distance_pc,
                aperture_au=aperture_au, label=label,
            )
        except Exception as e:
            print(f"  {name}: FAILED - {e}")
            continue

        row = sink_lookup[sink_id]
        L_acc = sinks.data[row, acclum_col]
        L_int = sinks.data[row, intlum_col]
        L_bol_true = L_acc + L_int

        sink_ids.append(sink_id)
        T_bols.append(bol["T_bol"])
        L_bol_obs_list.append(bol["L_bol"])
        L_acc_list.append(L_acc)
        L_int_list.append(L_int)
        L_bol_true_list.append(L_bol_true)

    L_bol_true_arr = np.array(L_bol_true_list)
    L_bol_obs_arr = np.array(L_bol_obs_list)

    table = {
        "sink_id": np.array(sink_ids),
        "T_bol": np.array(T_bols),
        "L_bol_obs": L_bol_obs_arr,
        "L_bol_true": L_bol_true_arr,
        "L_acc": np.array(L_acc_list),
        "L_int": np.array(L_int_list),
        "L_bol_ratio": np.where(L_bol_true_arr > 0,
                                 L_bol_obs_arr / L_bol_true_arr, 0.0),
    }

    print(f"\n  Processed {len(sink_ids)} sinks.")

    return table


# ============================================================
#  Spectral index and optically-thin mass
# ============================================================

def _planck_Bnu(freq_Hz, T_K):
    """Planck function B_ν(T) in W/m²/Hz/sr."""
    h = 6.62607015e-34
    k = 1.380649e-23
    c = 2.99792458e8
    x = h * freq_Hz / (k * T_K)
    return 2 * h * freq_Hz**3 / c**2 / (np.exp(x) - 1.0)


def compute_alpha_mass(
    pipeline_output,
    ramses_path=None,
    view="xy",
    distance_pc=140.0,
    wavelengths_mm=None,
    aperture_au=None,
    T_dust_K=20.0,
    kappa_ref_cm2g=2.3,
    kappa_ref_mm=1.3,
    kappa_beta=1.7,
    gas_to_dust_ratio=100.0,
    label="whole",
    sink_folders=None,
):
    """
    Compute spectral index α and optically-thin dust mass for each sink
    from multi-wavelength POLARIS images.

    α is fitted as a power law: log(F_ν) = α × log(ν) + const

    M_thin = F_ν × d² / (κ_ν × B_ν(T_dust))
    where κ_ν = kappa_ref × (ν / ν_ref)^beta

    Also computes M_true_aperture: the true (gas+dust) mass within
    the same aperture, integrated from the RADMC-3D density subbox.

    Parameters
    ----------
    pipeline_output : str or Path
        Pipeline output directory (e.g., test1/).
    ramses_path : str or Path or None
        RAMSES output directory. If provided, reads true sink masses.
    view : str
        Viewing direction: "xy", "xz", or "yz".
    distance_pc : float
        Source distance in parsecs.
    wavelengths_mm : list of float or None
        Wavelengths to use for α computation. If None, uses [0.87, 1.3, 3.0].
    aperture_au : float or None
        Aperture radius in AU. If None, integrate over the full image.
    T_dust_K : float
        Assumed dust temperature for optically-thin mass [K].
    kappa_ref_cm2g : float
        Dust opacity at the reference wavelength [cm²/g].
        Default: 2.3 cm²/g at 1.3 mm (Beckwith+ 1990, includes gas-to-dust 100).
    kappa_ref_mm : float
        Reference wavelength for kappa_ref [mm]. Default: 1.3 mm.
    kappa_beta : float
        Opacity power-law index: κ_ν ∝ ν^beta. Default: 1.7.
    gas_to_dust_ratio : float
        Gas-to-dust mass ratio for computing true total mass from the
        dust density grid. Default: 100.
    label : str
        Subdirectory label (default "whole").
    sink_folders : list of str or None
        If provided, only process these sinks.

    Returns
    -------
    table : dict
        Arrays keyed by: "sink_id", "alpha", "M_thin_Msun" (at each
        wavelength), "M_true_Msun", "fluxes_mJy", "wavelengths_mm".
    """
    from astropy.io import fits as pyfits

    if wavelengths_mm is None:
        wavelengths_mm = [0.87, 1.3, 3.0]

    pipeline_output = Path(pipeline_output)
    image_base = pipeline_output / "images" / "subboxes"

    c_light = 2.99792458e8
    pc2m_local = 3.085677581e16
    M_sun_kg = 1.989e30
    au2m_local = 1.4959787070e11

    distance_m = distance_pc * pc2m_local

    # Frequencies for the requested wavelengths
    freqs_Hz = np.array([c_light / (wl * 1e-3) for wl in wavelengths_mm])

    # Opacity at each wavelength: κ_ν = κ_ref × (ν / ν_ref)^β
    freq_ref = c_light / (kappa_ref_mm * 1e-3)
    kappa_cm2g = np.array([
        kappa_ref_cm2g * (f / freq_ref)**kappa_beta for f in freqs_Hz
    ])
    # Convert to SI: cm²/g → m²/kg (×0.1)
    kappa_SI = kappa_cm2g * 0.1

    # Find sink folders
    if sink_folders is not None:
        folder_names = sink_folders
    else:
        folder_names = sorted([
            d.name for d in image_base.iterdir()
            if d.is_dir() and d.name.startswith("sink_")
        ])

    # Read true masses from RAMSES if available
    sink_lookup = {}
    if ramses_path is not None:
        from radprocess.ramses import read as ramses_read
        sinks = ramses_read.sink_info(str(ramses_path))
        cols = sinks.columns
        id_col = cols.index("Id")
        m_col = cols.index("M[Msol]")
        for row_idx in range(sinks.num_sinks):
            sid = int(sinks.data[row_idx, id_col])
            sink_lookup[sid] = sinks.data[row_idx, m_col]

    # Process each sink
    sink_ids = []
    alphas = []
    all_fluxes_mJy = []
    M_thin_per_wl = []  # one mass estimate per wavelength
    M_sink_list = []
    M_true_aperture_list = []

    for name in folder_names:
        sink_id = int(name.replace("sink_", ""))
        data_dir = image_base / name / label / view / "data"

        if not data_dir.exists():
            print(f"  {name}: no {view} data, skipping")
            continue

        # Find FITS files and match to requested wavelengths
        all_fits = sorted(data_dir.glob("polaris_detector_nr*.fits*"))
        all_fits = [f for f in all_fits if "_sed" not in f.name]

        # Read wavelengths from all FITS headers
        fits_wavelengths = {}
        for ff in all_fits:
            hdul = pyfits.open(ff)
            h = hdul[0].header
            for key in ["WAVELENGTH1", "WAVELENGTH"]:
                if key in h:
                    wl_m = float(h[key])
                    fits_wavelengths[wl_m * 1e3] = ff  # key in mm
                    break
            hdul.close()

        # Match requested wavelengths to FITS files (within 1% tolerance)
        matched_files = []
        for wl_mm in wavelengths_mm:
            found = False
            for fits_wl_mm, ff in fits_wavelengths.items():
                if abs(fits_wl_mm - wl_mm) / wl_mm < 0.01:
                    matched_files.append(ff)
                    found = True
                    break
            if not found:
                matched_files.append(None)

        if any(f is None for f in matched_files):
            missing = [wl for wl, f in zip(wavelengths_mm, matched_files) if f is None]
            print(f"  {name}: missing wavelengths {missing} mm, skipping")
            continue

        # Extract fluxes
        fluxes_Jy = []
        for fits_file in matched_files:
            hdul = pyfits.open(fits_file)
            stokes_I = hdul[0].data[0, 0, :, :]  # Jy/px
            header = hdul[0].header
            hdul.close()

            cdelt_m = header["CDELT1"]
            ny, nx = stokes_I.shape

            if aperture_au is not None:
                cy_pix, cx_pix = ny / 2.0, nx / 2.0
                yy, xx = np.mgrid[0:ny, 0:nx]
                r_au = np.sqrt(
                    ((xx - cx_pix) * cdelt_m)**2 +
                    ((yy - cy_pix) * cdelt_m)**2
                ) / au2m_local
                mask = r_au <= aperture_au
                fluxes_Jy.append(np.sum(stokes_I[mask]))
            else:
                fluxes_Jy.append(np.sum(stokes_I))

        fluxes_Jy = np.array(fluxes_Jy)

        # Fit spectral index: log(F_ν) = α × log(ν) + const
        log_nu = np.log10(freqs_Hz)
        log_F = np.log10(np.maximum(fluxes_Jy, 1e-99))
        coeffs = np.polyfit(log_nu, log_F, 1)
        alpha = coeffs[0]

        # Optically-thin mass at each wavelength
        # M = F_ν × d² / (κ_ν × B_ν(T))
        fluxes_SI = fluxes_Jy * 1e-26  # Jy → W/m²/Hz
        M_thin = np.array([
            fluxes_SI[j] * distance_m**2 / (kappa_SI[j] * _planck_Bnu(freqs_Hz[j], T_dust_K))
            for j in range(len(wavelengths_mm))
        ])
        M_thin_Msun = M_thin / M_sun_kg

        # True mass from RAMSES sink particle
        M_sink = sink_lookup.get(sink_id, np.nan)

        # True mass within aperture from density subbox
        M_true_aperture = np.nan
        density_file = (pipeline_output / "radmc3d" / "subboxes" / name
                        / "dust_density_subbox.out")
        if density_file.exists():
            try:
                density_cube, bounds = _read_dust_density_subbox(density_file)
                xmin_d, xmax_d, ymin_d, ymax_d, zmin_d, zmax_d = bounds
                nx_d, ny_d, nz_d = density_cube.shape
                dx_cm = (xmax_d - xmin_d) / nx_d
                dy_cm = (ymax_d - ymin_d) / ny_d
                dz_cm = (zmax_d - zmin_d) / nz_d

                # Project along LOS → column density [g/cm²]
                # The density subbox is centered on the sink
                axis_for_proj = {"xy": "z", "xz": "y", "yz": "x"}[view]
                dlos = {"xy": dz_cm, "xz": dy_cm, "yz": dx_cm}[view]
                col_density = _project_density(density_cube, axis=axis_for_proj) * dlos

                # Pixel area in the image plane [cm²]
                if view == "xy":
                    pixel_area = dx_cm * dy_cm
                    img_nx, img_ny = nx_d, ny_d
                    half_x = (xmax_d - xmin_d) / 2.0
                    half_y = (ymax_d - ymin_d) / 2.0
                elif view == "xz":
                    pixel_area = dx_cm * dz_cm
                    img_nx, img_ny = nx_d, nz_d
                    half_x = (xmax_d - xmin_d) / 2.0
                    half_y = (zmax_d - zmin_d) / 2.0
                elif view == "yz":
                    pixel_area = dy_cm * dz_cm
                    img_nx, img_ny = ny_d, nz_d
                    half_x = (ymax_d - ymin_d) / 2.0
                    half_y = (zmax_d - zmin_d) / 2.0

                if aperture_au is not None:
                    # Aperture mask centered on the image center (= sink)
                    au2cm_local = 1.495978707e13
                    aperture_cm = aperture_au * au2cm_local
                    cx_pix = img_nx / 2.0
                    cy_pix = img_ny / 2.0
                    pix_dx = 2.0 * half_x / img_nx
                    pix_dy = 2.0 * half_y / img_ny
                    yy, xx = np.mgrid[0:col_density.shape[0], 0:col_density.shape[1]]
                    r_cm = np.sqrt(
                        ((xx - cx_pix) * pix_dx)**2 +
                        ((yy - cy_pix) * pix_dy)**2
                    )
                    mask_dens = r_cm <= aperture_cm
                    M_dust_g = np.sum(col_density[mask_dens]) * pixel_area
                else:
                    M_dust_g = np.sum(col_density) * pixel_area

                # Total mass = dust × (1 + gas-to-dust ratio)
                M_total_g = M_dust_g * (1.0 + gas_to_dust_ratio)
                M_sun_g = 1.989e33
                M_true_aperture = M_total_g / M_sun_g
            except Exception as e:
                print(f"    WARNING: could not compute true mass: {e}")

        sink_ids.append(sink_id)
        alphas.append(alpha)
        all_fluxes_mJy.append(fluxes_Jy * 1e3)  # Jy → mJy
        M_thin_per_wl.append(M_thin_Msun)
        M_sink_list.append(M_sink)
        M_true_aperture_list.append(M_true_aperture)

        print(f"  {name}: α = {alpha:.2f}, "
              f"F = [{', '.join(f'{f:.2f}' for f in fluxes_Jy*1e3)}] mJy, "
              f"M_thin(1.3mm) = {M_thin_Msun[1]:.4f} M_sun"
              f"{f', M_true(aperture) = {M_true_aperture:.4f} M_sun' if not np.isnan(M_true_aperture) else ''}")

    L_bol_true_arr = np.array(M_sink_list)
    M_true_ap_arr = np.array(M_true_aperture_list)
    M_thin_arr = np.array(M_thin_per_wl)

    table = {
        "sink_id": np.array(sink_ids),
        "alpha": np.array(alphas),
        "fluxes_mJy": np.array(all_fluxes_mJy),
        "wavelengths_mm": np.array(wavelengths_mm),
        "M_thin_Msun": M_thin_arr,
        "M_sink_Msun": L_bol_true_arr,
        "M_true_aperture_Msun": M_true_ap_arr,
        "M_ratio": np.where(M_true_ap_arr > 0,
                             M_thin_arr[:, 1] / M_true_ap_arr, np.nan),
        "T_dust_assumed_K": T_dust_K,
        "kappa_ref_cm2g": kappa_ref_cm2g,
        "kappa_beta": kappa_beta,
        "gas_to_dust_ratio": gas_to_dust_ratio,
        "aperture_au": aperture_au if aperture_au is not None else "full",
    }

    print(f"\n  Processed {len(sink_ids)} sinks.")
    print(f"  Assumed T_dust = {T_dust_K} K, "
          f"kappa = {kappa_ref_cm2g} cm²/g at {kappa_ref_mm} mm, beta = {kappa_beta}")
    print(f"  Gas-to-dust ratio = {gas_to_dust_ratio}")

    return table


# ============================================================
#  Mean optical depth within aperture
# ============================================================

def compute_tau_aperture(
    pipeline_output,
    view="xy",
    wavelength_idx=0,
    aperture_au=None,
    label="whole",
    sink_folders=None,
):
    """
    Compute the mean optical depth τ within an aperture for each sink.

    τ is read from the POLARIS FITS files (plane index 4).

    Parameters
    ----------
    pipeline_output : str or Path
        Pipeline output directory (e.g., test1/).
    view : str
        Viewing direction: "xy", "xz", or "yz".
    wavelength_idx : int
        Which wavelength (0-indexed).
    aperture_au : float or None
        Aperture radius in AU. If None, use the full image.
    label : str
        Subdirectory label (default "whole").
    sink_folders : list of str or None
        If provided, only process these sinks.

    Returns
    -------
    result : dict
        "sink_id": array of sink IDs
        "tau_mean": mean τ within aperture for each sink
        "tau_max": maximum τ within aperture
        "tau_median": median τ within aperture
    """
    from astropy.io import fits as pyfits

    pipeline_output = Path(pipeline_output)
    image_base = pipeline_output / "images" / "subboxes"
    au2m_local = 1.4959787070e11

    if sink_folders is not None:
        folder_names = sink_folders
    else:
        folder_names = sorted([
            d.name for d in image_base.iterdir()
            if d.is_dir() and d.name.startswith("sink_")
        ])

    sink_ids = []
    tau_means = []
    tau_maxs = []
    tau_medians = []
    tau_flux_weighteds = []
    skipped = 0

    for name in folder_names:
        sink_id = int(name.replace("sink_", ""))
        data_dir = image_base / name / label / view / "data"

        if not data_dir.exists():
            skipped += 1
            continue

        # Check how many wavelengths are available
        all_fits = sorted([
            f for f in data_dir.glob("polaris_detector_nr*.fits*")
            if "_sed" not in f.name
        ])

        if not all_fits:
            skipped += 1
            continue

        if wavelength_idx >= len(all_fits):
            if skipped == 0:  # print warning only once
                print(f"  WARNING: wavelength_idx={wavelength_idx} but only "
                      f"{len(all_fits)} wavelengths available (0-{len(all_fits)-1})")
            skipped += 1
            continue

        fits_file = all_fits[wavelength_idx]

        hdul = pyfits.open(fits_file)
        tau_map = hdul[0].data[4, 0, :, :]    # plane 4 = optical depth
        stokes_I = hdul[0].data[0, 0, :, :]   # plane 0 = Stokes I (Jy/px)
        header = hdul[0].header
        hdul.close()

        ny, nx = tau_map.shape
        cdelt_m = header["CDELT1"]

        if aperture_au is not None:
            cx_pix, cy_pix = nx / 2.0, ny / 2.0
            yy, xx = np.mgrid[0:ny, 0:nx]
            r_au = np.sqrt(
                ((xx - cx_pix) * cdelt_m)**2 +
                ((yy - cy_pix) * cdelt_m)**2
            ) / au2m_local
            mask = r_au <= aperture_au
            tau_in_ap = tau_map[mask]
            flux_in_ap = stokes_I[mask]
        else:
            tau_in_ap = tau_map.ravel()
            flux_in_ap = stokes_I.ravel()

        # Flux-weighted mean: τ_w = Σ(τ × I) / Σ(I)
        flux_sum = np.sum(flux_in_ap)
        if flux_sum > 0:
            tau_flux_weighted = np.sum(tau_in_ap * flux_in_ap) / flux_sum
        else:
            tau_flux_weighted = 0.0

        sink_ids.append(sink_id)
        tau_means.append(np.mean(tau_in_ap))
        tau_maxs.append(np.max(tau_in_ap))
        tau_medians.append(np.median(tau_in_ap))
        tau_flux_weighteds.append(tau_flux_weighted)

    if len(sink_ids) == 0:
        print(f"  WARNING: no sinks processed! Check wavelength_idx={wavelength_idx}. "
              f"Skipped {skipped} sinks.")

    # Read wavelength from the last opened file for the summary
    wl_label = ""
    if len(sink_ids) > 0:
        hdul = pyfits.open(fits_file)
        for key in ["WAVELENGTH1", "WAVELENGTH"]:
            if key in hdul[0].header:
                wl_m = float(hdul[0].header[key])
                wl_label = f", λ={wl_m*1e3:.2f} mm"
                break
        hdul.close()

    result = {
        "sink_id": np.array(sink_ids),
        "tau_mean": np.array(tau_means),
        "tau_max": np.array(tau_maxs),
        "tau_median": np.array(tau_medians),
        "tau_flux_weighted": np.array(tau_flux_weighteds),
    }

    print(f"  Computed τ for {len(sink_ids)} sinks "
          f"(view={view}, idx={wavelength_idx}{wl_label}, "
          f"aperture={aperture_au} AU)")

    return result