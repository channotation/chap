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
	dpi = args.dpi)

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
	dpi = args.dpi)

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
	dpi = args.dpi)

pl.close("pf_hydrophobicity_profile")

