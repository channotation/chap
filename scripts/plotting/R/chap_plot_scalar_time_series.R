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

# plot output size:
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

# colour:
line.colour <- "#737373"


################################################################################
# DATA READ-IN
################################################################################

# load first line from JSON file:
dat <- fromJSON(readLines(opt$filename, n = 1), flatten = FALSE)


################################################################################
# TIME SERIES PLOTS
################################################################################

# smoother length:
t.smooth <- 5000 # time scale of smoother
span.smooth <- t.smooth/abs(diff(range(as.data.frame(dat$pathwayScalarTimeSeries)$t)))


# Length
#-------------------------------------------------------------------------------

plt.length <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                     aes(x = t, 
                         y = length)) +
  geom_line(colour = line.colour) +
  geom_smooth(method = "loess",
              se = FALSE,
              colour = "gold",
              span = span.smooth,
              n = 2048,
              size = 1.2) +
  scale_x_continuous(expand = c(0,0)) + 
  xlab(expression(paste(t~~bgroup("(", ps,")")))) +
  ylab(expression(paste(L~~bgroup("(", nm,")")))) +
  ggtitle("Permeation Pathway Length over Time") +
  theme_chap

ggsave("time_series_pathway_length.png", 
       plt.length,
       width = plot.width,
       height = plot.height,
       dpi = opt$dpi)


# Volume
#-------------------------------------------------------------------------------

plt.volume <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                     aes(x = t, 
                         y = volume)) +
  geom_line(colour = line.colour) +
  geom_smooth(method = "loess",
              se = FALSE,
              colour = "gold",
              span = span.smooth,
              n = 2048,
              size = 1.2) +
  scale_x_continuous(expand = c(0,0)) + 
  xlab(expression(paste(t~~bgroup("(", ps,")")))) +
  ylab(expression(paste(V~~bgroup("(", nm^{3},")")))) +
  ggtitle("Permeation Pathway Volume over Time") +
  theme_chap

ggsave("time_series_pathway_volume.png", 
       plt.volume,
       width = plot.width,
       height = plot.height,
       dpi = opt$dpi)


# Number of Solvent Particles in Pathway
#-------------------------------------------------------------------------------

plt.solvent.number <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                             aes(x = t, 
                                 y = numPathway)) +
  geom_line(colour = line.colour) +
  geom_smooth(method = "loess",
              se = FALSE,
              colour = "gold",
              span = span.smooth,
              n = 2048,
              size = 1.2) +
  scale_x_continuous(expand = c(0,0)) + 
  xlab(expression(paste(t~~bgroup("(", ps,")")))) +
  ylab(expression(paste(N, phantom(bgroup("(",nm^{-3},")"))))) +
  ggtitle("Number of Solvent Particles in Pathway over Time") +
  theme_chap

ggsave("time_series_solvent_number_in_pathway.png", 
       plt.solvent.number,
       width = plot.width,
       height = plot.height,
       dpi = opt$dpi)


# Average Solvent Number Density in Pathway
#-------------------------------------------------------------------------------

plt.avg.solvent.number.density <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                                         aes(x = t, 
                                             y = numPathway/volume)) +
  geom_line(colour = line.colour) +
  geom_abline(intercept = 33.3679, 
              slope = 0, 
              linetype = 2, 
              colour = "darkblue",
              size = 0.75) + # bulk water density from Wikipedia
  geom_smooth(method = "loess",
              se = FALSE,
              colour = "gold",
              span = span.smooth,
              n = 2048,
              size = 1.2) +
  scale_x_continuous(expand = c(0,0)) + 
  xlab(expression(paste(t~~bgroup("(", ps,")")))) +
  ylab(expression(paste(N/V~~bgroup("(", nm^{-3},")")))) +
  ggtitle("Average Density in Pathway over Time") +
  theme_chap

ggsave("time_series_avg_solvent_number_density.png", 
       plt.avg.solvent.number.density,
       width = plot.width,
       height = plot.height,
       dpi = opt$dpi)


# Minimal Pore Radius
#-------------------------------------------------------------------------------

plt.min.pore.radius <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                              aes(x = t, 
                                  y = minRadius)) +
  geom_line(colour = line.colour) +
  geom_smooth(method = "loess",
              se = FALSE,
              colour = "gold",
              span = span.smooth,
              n = 2048,
              size = 1.2) +
  scale_x_continuous(expand = c(0,0)) + 
  xlab(expression(paste(t~~bgroup("(", ps,")")))) +
  ylab(expression(paste(min(R, s)~~bgroup("(", nm,")")))) +
  ggtitle("Minimal Pore Radius over Time") +
  theme_chap

ggsave("time_series_min_pore_radius.png", 
       plt.min.pore.radius,
       width = plot.width,
       height = plot.height,
       dpi = opt$dpi)


# Minimal Solvent Number Density
#-------------------------------------------------------------------------------

plt.min.solvent.density <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                                  aes(x = t, 
                                      y = minSolventDensity)) +
  geom_line(colour = line.colour) +
  geom_smooth(method = "loess",
              se = FALSE,
              colour = "gold",
              span = span.smooth,
              n = 2048,
              size = 1.2) +
  scale_x_continuous(expand = c(0,0)) + 
  xlab(expression(paste(t~~bgroup("(", ps,")")))) +
  ylab(expression(paste(min(n, s)~~bgroup("(", nm^{-3},")")))) +
  ggtitle("Minimal Solvent Number Density over Time") +
  theme_chap

ggsave("time_series_min_solvent_number_density.png", 
       plt.min.solvent.density,
       width = plot.width,
       height = plot.height,
       dpi = opt$dpi)


# Location of Minimal Pore Radius
#-------------------------------------------------------------------------------

plt.argmin.pore.radius <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                                 aes(x = t, 
                                     y = argMinRadius)) +
  geom_point(colour = line.colour) +
  scale_x_continuous(expand = c(0,0)) + 
  xlab(expression(paste(t~~bgroup("(", ps,")")))) +
  ylab(expression(paste(arg*min(R, s)~~bgroup("(", nm,")")))) +
  ggtitle("Location of Minimal Pore Radius over Time") +
  theme_chap

ggsave("time_series_argmin_pore_radius.png", 
       plt.argmin.pore.radius,
       width = plot.width,
       height = plot.height,
       dpi = opt$dpi)


# Location of Minimal Solvent Number Density
#-------------------------------------------------------------------------------

plt.argmin.solvent.density <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                                     aes(x = t, 
                                         y = argMinSolventDensity)) +
  geom_point(colour = line.colour) +
  scale_x_continuous(expand = c(0,0)) + 
  xlab(expression(paste(t~~bgroup("(", ps,")")))) +
  ylab(expression(paste(arg*min(n, s)~~bgroup("(", nm^{-3},")")))) +
  ggtitle("Location of Minimal Density over Time") +
  theme_chap

ggsave("time_series_argmin_solvent_number_density.png", 
       plt.argmin.solvent.density,
       width = plot.width,
       height = plot.height,
       dpi = opt$dpi)
