#!/usr/bin/env python2

# CHAP - The Channel Annotation Package
# 
# Copyright (c) 2016 - 2018 Gianni Klesse, Shanlin Rao, Mark S. P. Sansom, and 
# Stephen J. Tucker
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


################################################################################
# CONFIGURATION
################################################################################

# load libraries:
import json                             # read in JSON files
import numpy as np                      # manipulate numeric vectors
from matplotlib import pyplot as pl     # plotting facilities
import argparse                         # parse command line arguments

# get parameters from user input:
parser = argparse.ArgumentParser()
parser.add_argument(
    "-filename",
    nargs = "?",
    const = "output.json",
    default = "output.json")
parser.add_argument("-dpi",
    nargs = "?",
    const = 1200,
    default = 1200,
    type = int)
args = parser.parse_args()


################################################################################
# DATA READ-IN
################################################################################

# load output data from JSON file:
with open(args.filename) as data_file:
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
	dpi = args.dpi)

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
	dpi = args.dpi)

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
pl.ylabel("n (nm$^{-3}$)")

pl.savefig(
    "time_averaged_solvent_number_density_profile.png",
	dpi = args.dpi)

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
	dpi = args.dpi)

pl.close("energy_profile")

