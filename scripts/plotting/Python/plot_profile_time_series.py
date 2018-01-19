#!/usr/bin/env python2

################################################################################
# SETTINGS (CHANGE FILENAME HERE)
################################################################################

# name of data file:
filename = "output.json"

# plot output parameters:
plot_dpi = 300


################################################################################
# CONFIGURATION
################################################################################

# load libraries:
import json
import numpy as np
from matplotlib import pyplot as pl


################################################################################
# DATA READ-IN
################################################################################

# load output data from JSON file:
with open(filename) as data_file:
    data = json.load(data_file)


################################################################################
# PATHWAY PROFILE PLOTS
################################################################################


# Radius Profile Time Series
#-------------------------------------------------------------------------------

pl.figure("radius_profile")

s = np.array(data["pathwayProfileTimeSeries"]["s"])
t = np.array(data["pathwayProfileTimeSeries"]["t"])
r = np.array(data["pathwayProfileTimeSeries"]["radius"])

num_t = np.size(np.unique(t))
num_s = np.size(np.unique(s))

S = s.reshape(num_t, num_s)
T = t.reshape(num_t, num_s)
R = r.reshape(num_t, num_s)

pl.pcolormesh(
    T,
    S,
    R,
    cmap = "YlOrBr_r")

cbar = pl.colorbar()
cbar.ax.set_ylabel("R (nm)")

pl.title("Radius Profile over Time")
pl.xlabel("t (ps)")
pl.ylabel("s (nm)")

pl.savefig(
    "time_series_radius_profile.png",
	dpi = plot_dpi)

pl.close("radius_profile")


# Solvent Density Profile Time Series
#-------------------------------------------------------------------------------

pl.figure("density_profile")

s = np.array(data["pathwayProfileTimeSeries"]["s"])
t = np.array(data["pathwayProfileTimeSeries"]["t"])
n = np.array(data["pathwayProfileTimeSeries"]["density"])

num_t = np.size(np.unique(t))
num_s = np.size(np.unique(s))

S = s.reshape(num_t, num_s)
T = t.reshape(num_t, num_s)
N = n.reshape(num_t, num_s)

pl.pcolormesh(
    T,
    S,
    N,
    cmap = "Blues")

cbar = pl.colorbar()
cbar.ax.set_ylabel("n ($\mathrm{nm}^{-3}$)")

pl.title("Number Density Profile over Time")
pl.xlabel("t (ps)")
pl.ylabel("s (nm)")

pl.savefig(
    "time_series_number_density_profile.png",
	dpi = plot_dpi)

pl.close("density_profile")


# Pore-facing Hydrophobicity Profile Time Series
#-------------------------------------------------------------------------------

pl.figure("pf_hydrophobicity_profile")

s = np.array(data["pathwayProfileTimeSeries"]["s"])
t = np.array(data["pathwayProfileTimeSeries"]["t"])
h = np.array(data["pathwayProfileTimeSeries"]["pfHydrophobicity"])

num_t = np.size(np.unique(t))
num_s = np.size(np.unique(s))

S = s.reshape(num_t, num_s)
T = t.reshape(num_t, num_s)
H = h.reshape(num_t, num_s)

pl.pcolormesh(
    T,
    S,
    H,
    cmap = "BrBG_r")
pl.clim(
    -max(abs(np.array(data["pathwayProfileTimeSeries"]["pfHydrophobicity"]))),
    max(abs(np.array(data["pathwayProfileTimeSeries"]["pfHydrophobicity"]))))

cbar = pl.colorbar()
cbar.ax.set_ylabel("H (a.u.)")

pl.title("Hydrophobicity Profile over Time")
pl.xlabel("t (ps)")
pl.ylabel("s (nm)")

pl.savefig(
    "time_series_hydrophobicity_profile.png",
	dpi = plot_dpi)

pl.close("pf_hydrophobicity_profile")

