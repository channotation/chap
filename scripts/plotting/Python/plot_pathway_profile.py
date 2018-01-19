#!/usr/bin/env python2

################################################################################
# SETTINGS (CHANGE FILENAME HERE)
################################################################################

# name of data file:
filename = "output.json"

# plot output parameters:
plot_dpi = 1200


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


# Radius Profile with Residue Positions
#-------------------------------------------------------------------------------

pl.figure("radius_profile")

pl.plot(
	np.array(data["pathwayProfile"]["s"]),
    np.array(data["pathwayProfile"]["radiusMean"]),
	"k-")

pl.fill_between(
	np.array(data["pathwayProfile"]["s"]),
	np.array(data["pathwayProfile"]["radiusMin"]),
	np.array(data["pathwayProfile"]["radiusMax"]),
	facecolor = "#000000",
	alpha = 0.1)

radius_sd = np.array(data["pathwayProfile"]["radiusSd"])
pl.fill_between(
	np.array(data["pathwayProfile"]["s"]),
	np.array(data["pathwayProfile"]["radiusMean"]) - radius_sd,
    np.array(data["pathwayProfile"]["radiusMean"]) + radius_sd,
	facecolor = "#000000",
	alpha = 0.2)

pf = np.array(data["residueSummary"]["poreFacing"]["mean"]) > 0.5
pl.scatter(
	np.array(data["residueSummary"]["s"]["mean"])[pf],
	np.array(data["residueSummary"]["rho"]["mean"])[pf],
    c =	np.array(data["residueSummary"]["hydrophobicity"])[pf],
    marker = "o",
    cmap = "BrBG_r")
pl.clim(
    -max(abs(np.array(data["residueSummary"]["hydrophobicity"]))),
    max(abs(np.array(data["residueSummary"]["hydrophobicity"]))))
cbar = pl.colorbar()
cbar.ax.set_ylabel("Hydrophobicity (a.u.)")

pl.margins(x = 0)
pl.title("Time-Averaged Radius Profile")
pl.xlabel("s (nm)")
pl.ylabel("R (nm)")

pl.savefig(
    "time_averaged_radius_profile.png",
	dpi = plot_dpi)

pl.close("radius_profile")


# Hydrophobicity Profile with Residue Positions
#-------------------------------------------------------------------------------

pl.figure("hydrophobicity_profile")

pl.plot(
	np.array(data["pathwayProfile"]["s"]),
    np.array(data["pathwayProfile"]["pfHydrophobicityMean"]),
	"k-")

pl.fill_between(
	np.array(data["pathwayProfile"]["s"]),
	np.array(data["pathwayProfile"]["pfHydrophobicityMin"]),
	np.array(data["pathwayProfile"]["pfHydrophobicityMax"]),
	facecolor = "#000000",
	alpha = 0.1)

hydrophobicity_sd = np.array(data["pathwayProfile"]["pfHydrophobicitySd"])
pl.fill_between(
	np.array(data["pathwayProfile"]["s"]),
	np.array(data["pathwayProfile"]["pfHydrophobicityMean"]) - hydrophobicity_sd,
    np.array(data["pathwayProfile"]["pfHydrophobicityMean"]) + hydrophobicity_sd,
	facecolor = "#000000",
	alpha = 0.2)

pf = np.array(data["residueSummary"]["poreFacing"]["mean"]) > 0.5
pl.scatter(
	np.array(data["residueSummary"]["s"]["mean"])[pf],
	np.array(data["residueSummary"]["hydrophobicity"])[pf],
    c =	np.array(data["residueSummary"]["hydrophobicity"])[pf],
    marker = "o",
    cmap = "BrBG_r")
pl.clim(
    -max(abs(np.array(data["residueSummary"]["hydrophobicity"]))),
    max(abs(np.array(data["residueSummary"]["hydrophobicity"]))))
cbar = pl.colorbar()
cbar.ax.set_ylabel("Hydrophobicity (a.u.)")

pl.margins(x = 0)
pl.title("Time-Averaged Hydrophobicity Profile")
pl.xlabel("s (nm)")
pl.ylabel("H (a.u.)")

pl.savefig(
    "time_averaged_hydrophobicity_profile.png",
	dpi = plot_dpi)

pl.close("hydrophobicity_profile")


# Solvent Number Density Profile
#-------------------------------------------------------------------------------

pl.figure("density_profile")

pl.axhline(
	y = 33.3679,
	linestyle = "dashed")

pl.plot(
	np.array(data["pathwayProfile"]["s"]),
    np.array(data["pathwayProfile"]["densityMean"]),
	"k-")

pl.fill_between(
	np.array(data["pathwayProfile"]["s"]),
	np.array(data["pathwayProfile"]["densityMin"]),
	np.array(data["pathwayProfile"]["densityMax"]),
	facecolor = "#000000",
	alpha = 0.1)

density_sd = np.array(data["pathwayProfile"]["densitySd"])
pl.fill_between(
	np.array(data["pathwayProfile"]["s"]),
	np.array(data["pathwayProfile"]["densityMean"]) - density_sd,
    np.array(data["pathwayProfile"]["densityMean"]) + density_sd,
	facecolor = "#000000",
	alpha = 0.2)

pl.margins(x = 0)
pl.title("Time-Averaged Solvent Number Density Profile")
pl.xlabel("s (nm)")
pl.ylabel("H (a.u.)")

pl.savefig(
    "time_averaged_solvent_number_density_profile.png",
	dpi = plot_dpi)

pl.close("density_profile")


# Solvent Number Density Profile
#-------------------------------------------------------------------------------

pl.figure("energy_profile")

pl.plot(
	np.array(data["pathwayProfile"]["s"]),
    np.array(data["pathwayProfile"]["energyMean"]),
	"k-")

energy_sd = np.array(data["pathwayProfile"]["energySd"])
pl.fill_between(
	np.array(data["pathwayProfile"]["s"]),
	np.array(data["pathwayProfile"]["energyMean"]) - energy_sd,
    np.array(data["pathwayProfile"]["energyMean"]) + energy_sd,
	facecolor = "#000000",
	alpha = 0.2)

pl.margins(x = 0)
pl.title("Time-Averaged Solvent Number Density Profile")
pl.xlabel("s (nm)")
pl.ylabel("G (kT)")

pl.savefig(
    "time_averaged_free_energy_profile.png",
	dpi = plot_dpi)

pl.close("energy_profile")

