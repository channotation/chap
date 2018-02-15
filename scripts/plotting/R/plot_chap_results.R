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
# SETTINGS
################################################################################

# libraries:
library(ggplot2)  # plotting
library(jsonlite) # parsing JSON files

# working directory:
setwd("~/repos/chap/bin")

# name of data file:
filename <- "output.json"


################################################################################
# DATA READ-IN
################################################################################

# load first line from JSON file:
dat <- fromJSON(readLines(filename, n = 1), flatten = FALSE)


################################################################################
# PLOTTING
################################################################################

# Radius Profile
#-------------------------------------------------------------------------------

plt.radius <- ggplot(as.data.frame(dat$pathProfile),
                     aes(x = s,
                         y = radiusMean,
                         ymin = radiusMean - radiusSd,
                         ymax = radiusMean + radiusSd)) +
  geom_ribbon(alpha = 0.2,
              aes(ymin = radiusMin,
                  ymax = radiusMax)) +
  geom_ribbon(alpha = 0.4) +
  geom_line() +
  xlab(expression(paste(s~~bgroup("(",nm,")")))) +
  ylab(expression(paste(R~~bgroup("(",nm,")")))) +
  ggtitle("Time Averaged Radius Profile")

plt.radius

# save plot to as files:
ggsave("radius_profile.png", plt.radius)
ggsave("radius_profile.pdf", plt.radius)


# Pore Hydration / Solvent Density Profile
#-------------------------------------------------------------------------------

plt.density <- ggplot(as.data.frame(dat$outfile$pathProfile),
                      aes(x = s,
                          y = densityMean,
                          ymin = densityMean - densitySd,
                          ymax = densityMean + densitySd)) +
  geom_ribbon(alpha = 0.2,
              aes(ymin = densityMin,
                  ymax = densityMax)) +
  geom_ribbon(alpha = 0.4) +
  geom_line() +
  xlab(expression(paste(s~~bgroup("(",nm,")")))) +
  ylab(expression(paste(n~~bgroup("(",nm^{-3},")")))) +
  ggtitle("Time Averaged Solvent Number Density")

plt.density

# save plot as file:
ggsave("solvent_number_density.png", plt.density)
ggsave("solvent_number_density.pdf", plt.density)


# Energy Profile
#-------------------------------------------------------------------------------

plt.energy <- ggplot(as.data.frame(dat$pathProfile),
                     aes(x = s,
                         y = energyMean)) +
  geom_ribbon(alpha = 0.4,
              aes(ymin = energyMean - energySd,
                  ymax = energyMean + energySd)) +
  geom_ribbon(alpha = 0.2,
              aes(ymin = energyMin,
                  ymax = energyMax)) +
  geom_line() +
  xlab(expression(paste(s~~bgroup("(",nm,")")))) +
  ylab(expression(paste(G~~bgroup("(",k[B]*T,")")))) +
  ggtitle("Solvent Free Energy Profile")

plt.energy

# save plot to file:
ggsave("energy_profile.png", plt.energy)
ggsave("energy_profile.pdf", plt.energy)
  








