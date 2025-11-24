#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt

from ..constants.constants import au2cm, M_sun, R_sun


def plot_sink_columns(sink_data, x_col, y_col):
    """
    sink_data: list of dicts from sink_info
    x_col, y_col: strings
    """
    x = [float(row[x_col]) for row in sink_data]
    y = [float(row[y_col]) for row in sink_data]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(x, y)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title(f"{y_col} vs {x_col}")
    return fig
