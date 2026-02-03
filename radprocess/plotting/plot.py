#!/usr/bin/env python
# -*- coding: utf-8 -*-
import math

import numpy as np
import matplotlib.pyplot as plt

from ..constants.constants import au2cm, M_sun, R_sun


def plot_sink_columns(sink_data, x_col, y_col):
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



def plot_ramses_histogram(root, y_field):
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
