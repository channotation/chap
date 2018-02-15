#!/usr/bin/env Rscript

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
if( !require(optparse) )
{
  install.packages("optparse", repos = "http://cran.us.r-project.org")
  library(optparse)
}
if( !require(jsonlite) )
{
  install.packages("jsonlite", repos = "http://cran.us.r-project.org")
  library(jsonlite)
}
if( !require(ggplot2) )
{
  install.packages("ggplot2", repos = "http://cran.us.r-project.org")
  library(ggplot2)
}

# get command line options:
option_list = list(
  make_option(c("--filename"),
              action = "store",
              default = "output.json",
              type = "character",
              help = "Name of the input file."),
  make_option(c("--dpi"),
              action = "store",
              default = 300,
              type = "numeric",
              help = "Resolution of plot in dots per inch.")
)
opt = parse_args(OptionParser(option_list=option_list))


################################################################################
# PLOT APPEARANCE
################################################################################

# plot output parameters:
plot.width.cm <- 5.8
plot.height.cm <- plot.width.cm * 0.75
plot.width <- unit(plot.width.cm, "cm")
plot.height <- unit(plot.height.cm, "cm")

# plot appearance:
theme_chap <- theme(text = element_text(size = 13),
                    axis.text.x = element_text(size = 12,
                                               colour = "black"),
                    axis.text.y = element_text(size = 12,
                                               colour = "black"),
                    panel.border = element_rect(colour = "black", 
                                                fill = NA, 
                                                size=1),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank()) 


################################################################################
# DATA READ-IN
################################################################################

# load first line from JSON file:
dat <- fromJSON(readLines(opt$filename, n = 1), flatten = FALSE)


################################################################################
# PATHWAY PROFILE PLOTS
################################################################################


# Radius Profile with Residue Positions
#-------------------------------------------------------------------------------

data <- as.data.frame(dat$residueSummary)
plt.radius.profile <- ggplot(data = data[dat$residueSummary$poreFacing$mean >= 0.1 &
                                         dat$residueSummary$id < length(dat$residueSummary$id),],) +
  geom_line(data = as.data.frame(dat$pathwayProfile), 
            aes(x = s, 
                y = radiusMean)) +
  geom_ribbon(data = as.data.frame(dat$pathwayProfile), 
              aes(x = s,
                  ymin = radiusMin, 
                  ymax = radiusMax), 
              alpha = 0.1) +
  geom_ribbon(data = as.data.frame(dat$pathwayProfile), 
              aes(x = s,
                  ymin = radiusMean - radiusSd, 
                  ymax = radiusMean + radiusSd), 
              alpha = 0.2) +
  geom_point(aes(x = s.mean,
                 y = rho.mean,
                 colour = hydrophobicity),
             size = 4) +
  scale_colour_distiller(palette = "BrBG",
                         name = expression(paste(H~~bgroup("(",a.u.,")"))),
                         limits = c(-max(abs(data$hydrophobicity)), 
                                    max(abs(data$hydrophobicity))),
                         guide = guide_colourbar(barheight = 10,
                                                 barwidth = 1.3)) +
  scale_x_continuous(expand = c(0, 0),
                     name = expression(paste(s~~bgroup("(",nm,")")))) +
  scale_y_continuous(name = expression(paste(R~~bgroup("(",nm,")")))) +
  ggtitle("Time-Averaged Radius Profile") +
  theme_chap

ggsave("time_averaged_radius_profile.png", 
       plt.radius.profile,
       width = plot.width,
       height = plot.height,
       dpi = opt$dpi)


# Hydrophobicity Profile with Residue Positions
#-------------------------------------------------------------------------------

data <- as.data.frame(dat$residueSummary)
plt.hydrophobicity <- ggplot(as.data.frame(dat$pathwayProfile),
                             aes(x = s,
                                 y = pfHydrophobicityMean)) +
  geom_line() +
  geom_point(data = data[data$poreFacing.mean > 0.5,],
             aes(x = s.mean,
                 y = hydrophobicity,
                 colour = hydrophobicity),
             size = 4) +
  geom_ribbon(aes(ymin = pfHydrophobicityMin, 
                  ymax = pfHydrophobicityMax), 
              alpha = 0.1) +
  geom_ribbon(aes(ymin = pfHydrophobicityMean - pfHydrophobicitySd, 
                  ymax = pfHydrophobicityMean + pfHydrophobicitySd), 
              alpha = 0.2) +
  scale_x_continuous(expand = c(0, 0),
                     name = expression(paste(s~~bgroup("(",nm,")")))) +
  scale_y_continuous(name = expression(paste(H~~bgroup("(",nm,")")))) +
  scale_colour_distiller(palette = "BrBG",
                         name = expression(paste(H~~bgroup("(",a.u.,")"))),
                         limits = c(-max(abs(data$hydrophobicity)), 
                                    max(abs(data$hydrophobicity))),
                         guide = guide_colourbar(barheight = 10,
                                                 barwidth = 1.3)) +
  ggtitle("Time-Averaged Hydrophobicity") +
  theme_chap

ggsave("time_averaged_hydrophobicity_profile.png", 
       plt.hydrophobicity,
       width = plot.width,
       height = plot.height,
       dpi = opt$dpi)


# Solvent Number Density Profile
#-------------------------------------------------------------------------------

plt.solvent.number.density <- ggplot(as.data.frame(dat$pathwayProfile),
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
  scale_x_continuous(expand = c(0, 0),
                     name = expression(paste(s~~bgroup("(",nm,")")))) +
  scale_y_continuous(name = expression(paste(n~~bgroup("(",nm^{-3},")")))) +
  ggtitle("Time-Averaged Solvent Number Density Profile") +
  theme_chap

ggsave("time_averaged_solvent_number_density_profile.png", 
       plt.solvent.number.density,
       width = plot.width,
       height = plot.height,
       dpi = opt$dpi)


# Solvent Free Energy Profile
#-------------------------------------------------------------------------------

plt.free.energy <- ggplot(as.data.frame(dat$pathwayProfile),
                          aes(x = s,
                              y = energyMean)) +
  geom_line() +
  geom_ribbon(aes(ymin = energyMean - energySd, 
                  ymax = energyMean + energySd), 
              alpha = 0.2) +
  scale_x_continuous(expand = c(0, 0),
                     name = expression(paste(s~~bgroup("(",nm,")")))) +
  scale_y_continuous(name = expression(paste(G~~bgroup("(",k[B]*T,")")))) +
  ggtitle("Time-Averaged Solvent Free Energy Profile") +
  theme_chap

ggsave("time_averaged_free_energy_profile.png", 
       plt.free.energy,
       width = plot.width,
       height = plot.height,
       dpi = opt$dpi)
