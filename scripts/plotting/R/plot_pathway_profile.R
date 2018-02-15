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
# SETTINGS (CHANGE FILENAME HERE)
################################################################################

# libraries:
library(ggplot2)   # plotting
library(ggrepel)   # improved labels
library(gridExtra) # plot panel arrangement
library(jsonlite)  # parsing JSON files

# name of data file:
filename <- "output.json"

# plot output parameters:
plot.width.cm <- 21.0/2
plot.height.cm <- 29.7/4

# plot appearance:
theme_chap <- theme_grey(base_size = 18)


################################################################################
# DATA READ-IN
################################################################################

# load first line from JSON file:
dat <- fromJSON(readLines(filename, n = 1), flatten = FALSE)


################################################################################
# PATHWAY PROFILE PLOTS
################################################################################

# individual plot size:
plot.width <- unit(plot.width.cm, "cm")
plot.height <- unit(plot.height.cm, "cm")


# Radius Profile
#-------------------------------------------------------------------------------

data <- as.data.frame(dat$residueSummary)
plt.radius.profile <- ggplot(data = data[dat$residueSummary$poreFacing$mean >= 0.1 &
                                         dat$residueSummary$id < length(dat$residueSummary$id),],) +
  geom_line(data = as.data.frame(dat$pathProfile), 
            aes(x = s, 
                y = radiusMean)) +
  geom_ribbon(data = as.data.frame(dat$pathProfile), 
              aes(x = s,
                  ymin = radiusMin, 
                  ymax = radiusMax), 
              alpha = 0.1) +
  geom_ribbon(data = as.data.frame(dat$pathProfile), 
              aes(x = s,
                  ymin = radiusMean - radiusSd, 
                  ymax = radiusMean + radiusSd), 
              alpha = 0.2) +
  geom_point(aes(x = s.mean,
                 y = rho.mean,
                 colour = hydrophobicity),
             size = 2) +
  scale_colour_distiller(palette = "RdBu",
                         name = expression(paste("hydrophobicity")),
                         limits = c(-max(abs(data$hydrophobicity)), 
                                    max(abs(data$hydrophobicity)))) +
  xlab(expression(paste(s~~bgroup("(",nm,")")))) +
  ylab(expression(paste(R~~bgroup("(",nm,")")))) +
  ggtitle("Time-Averaged Pore-Facing Residue Positions") +
  theme_chap

ggsave("time_averaged_radius_profile.pdf", 
       plt.radius.profile,
       width = plot.width,
       height = plot.height)


# Radius Profile with Residue Positions
#-------------------------------------------------------------------------------

data <- as.data.frame(dat$residueSummary)
plt.hydrophobicity <- ggplot(as.data.frame(dat$pathProfile),
                             aes(x = s,
                                 y = pfHydrophobicityMean)) +
  geom_line() +
  geom_point(data = data[data$poreFacing.mean > 0.1,],
             aes(x = s.mean,
                 y = hydrophobicity,
                 colour = name),
             size = 2) +
  geom_ribbon(aes(ymin = pfHydrophobicityMin, ymax = pfHydrophobicityMax), alpha = 0.1) +
  geom_ribbon(aes(ymin = pfHydrophobicityMean - pfHydrophobicitySd, ymax = pfHydrophobicityMean + pfHydrophobicitySd), alpha = 0.2) +
  ggtitle("Time-Averaged Hydrophobicity due to Pore-Facing Residues") +
  xlab(expression(paste(s~~bgroup("(",nm,")")))) +
  ylab(expression(paste("hydrophobicity"))) +
  theme_chap

ggsave("time_averaged_hydrophobicity.pdf", 
       plt.hydrophobicity,
       width = plot.width,
       height = plot.height)


# Solvent Number Density Profile
#-------------------------------------------------------------------------------

plt.solvent.number.density <- ggplot(as.data.frame(dat$pathProfile),
                                     aes(x = s,
                                         y = densityMean)) +
  geom_abline(intercept = 33.3679, 
              slope = 0, 
              linetype = 2, 
              colour = "darkblue") + # bulk water density from Wikipedia
  geom_line() +
  geom_ribbon(aes(ymin = densityMin, 
                  ymax = densityMax), 
              alpha = 0.1) +
  geom_ribbon(aes(ymin = densityMean - densitySd, 
                  ymax = densityMean + densitySd), 
              alpha = 0.2) +
  xlab(expression(paste(s~~bgroup("(",nm,")")))) +
  ylab(expression(paste(n~~bgroup("(",nm^{-3},")")))) +
  ggtitle("Time-Averaged Solvent Number Density Profile") +
  theme_chap

ggsave("time_averaged_solvent_number_density_profile.pdf", 
       plt.solvent.number.density,
       width = plot.width,
       height = plot.height)


# Solvent Free Energy Profile
#-------------------------------------------------------------------------------

plt.free.energy <- ggplot(as.data.frame(dat$pathProfile),
                          aes(x = s,
                              y = energyMean)) +
  geom_line() +
  geom_ribbon(aes(ymin = energyMin, 
                  ymax = energyMax), 
              alpha = 0.1) +
  geom_ribbon(aes(ymin = energyMean - energySd, 
                  ymax = energyMean + energySd), 
              alpha = 0.2) +
  xlab(expression(paste(s~~bgroup("(",nm,")")))) +
  ylab(expression(paste(G~~bgroup("(",kT[B],")")))) +
  ggtitle("Time-Averaged Solvent Solvent Free Energy Profile") +
  theme_chap

ggsave("time_averaged_free_energy_profile.pdf", 
       plt.free.energy,
       width = plot.width,
       height = plot.height)


# Combined Plot
#-------------------------------------------------------------------------------

plt.profile <- grid.arrange(plt.radius.profile, 
                            plt.hydrophobicity,
                            plt.solvent.number.density,
                            plt.free.energy,
                            ncol = 2)

ggsave("pathway_profiles.pdf", 
       plt.profile, 
       width = unit(2*plot.width.cm, "cm"), 
       height = unit(2*plot.height.cm, "cm"))

