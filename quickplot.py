#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import itertools


COLOR_CELLS = (34.0/255.0, 139.0/255.0, 34/255.0)
COLOR_SETPOINT = "red"
COLOR_FILL_CELLS = COLOR_CELLS
COLOR_FILL_DIFF = "orange"
ALPHA_FILL = 0.6

PLOT_MULTIPLIER_SIZE_X = 4.0
PLOT_MULTIPLIER_SIZE_Y = 3.0


# Read information from csv file
df = pd.read_csv("information.csv")

# Get unique ids of contorllers
id_x = sorted(np.unique(df.id_x))
id_y = sorted(np.unique(df.id_y))

# Determina all possible combinations
combinations = itertools.product(id_x, id_y)

# Create a figure which has more subplots than needed
fig, ax = plt.subplots(max(id_x)+1, max(id_y)+1, figsize=(max(id_x)*PLOT_MULTIPLIER_SIZE_X, max(id_y)*PLOT_MULTIPLIER_SIZE_Y), sharex=True, sharey=True)

# Plot all possible combinations
for i, j in combinations:
    # Get entries which match the indices
    entries = df[df.id_x == i][df.id_y == j]
    cells = entries.N_cells
    setpoint = entries.setpoint

    # Filter out if we have a dummy controller who is only used for measuring and not controlling
    if not sum(np.isnan(setpoint)):
        ax[i, j].plot(cells, label="Controller ({:1},{:1})".format(i, j), color=COLOR_CELLS)
        ax[i, j].fill_between(entries.index, cells, setpoint, color=COLOR_FILL_DIFF, alpha=ALPHA_FILL)
        ax[i, j].legend()
        ax[i, j].plot(setpoint, color=COLOR_SETPOINT)
    else:
        ax[i, j].plot(cells, label="Non-controlled voxel", color=COLOR_CELLS)
        ax[i, j].fill_between(entries.index, cells, 0, color=COLOR_FILL_CELLS, alpha=ALPHA_FILL)
        ax[i, j].legend()

# Adjust whitespaces between plots such that there are no gaps
fig.subplots_adjust(wspace=0, hspace=0)
fig.suptitle("Time evolution of Control Regions")
fig.supxlabel("Time [min]")
fig.supylabel("Number of cells")
fig.savefig("controller_regions_plots.svg")
plt.show()
