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
theme_heatmap <- theme(text = element_text(size = 13),
                       axis.text.x = element_text(size = 12,
                                                  colour = "black"),
                       axis.text.y = element_text(size = 12,
                                                  colour = "black"),
                       panel.border = element_rect(colour = "black", 
                                                   fill=NA, 
                                                   size=1)) 


################################################################################
# DATA READ-IN
################################################################################

# load first line from JSON file:
dat <- fromJSON(readLines(opt$filename, n = 1), flatten = FALSE)


################################################################################
# PATHWAY PROFILE PLOTS
################################################################################


# Radius Profile Time Series
#-------------------------------------------------------------------------------

plt.radius.profile.ts <- ggplot(data = as.data.frame(dat$pathwayProfileTimeSeries),
                                aes(x = t,
                                    y = s,
                                    fill = radius)) +
  geom_tile() +
  scale_fill_distiller(palette = "YlOrBr",
                       direction = -1,
                       name = expression(paste(R~~bgroup("(", nm,")"))),
                       guide = guide_colourbar(barheight = 10,
                                               barwidth = 1.3)) +
  scale_x_continuous(expand = c(0, 0),
                     name = expression(paste(t~~bgroup("(", ps,")")))) +
  scale_y_continuous(expand = c(0, 0),
                     name = expression(paste(s~~bgroup("(", nm,")")))) +
  ggtitle("Radius Profile over Time") +
  theme_heatmap

ggsave("time_series_radius_profile.png", 
       plt.radius.profile.ts,
       width = plot.width,
       height = plot.height,
       dpi = opt$dpi)


# Solvent Density Profile Time Series
#-------------------------------------------------------------------------------

plt.density.ts <- ggplot(data = as.data.frame(dat$pathwayProfileTimeSeries),
                         aes(x = t,
                             y = s,
                             fill = density)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues",
                       direction = 1,
                       name = expression(paste(n~~bgroup("(", nm^{-3},")"))),
                       guide = guide_colourbar(barheight = 10,
                                               barwidth = 1.3)) +
  scale_x_continuous(expand = c(0, 0),
                     name = expression(paste(t~~bgroup("(", ps,")")))) +
  scale_y_continuous(expand = c(0, 0),
                     name = expression(paste(s~~bgroup("(", nm,")")))) +
  ggtitle("Number Density Profile over Time") +
  theme_heatmap
  
ggsave("time_series_number_density_profile.png", 
       plt.density.ts,
       width = plot.width,
       height = plot.height,
       dpi = opt$dpi)


# Pore-facing Hydrophobicity Profile Time Series
#-------------------------------------------------------------------------------

plt.pf.hydrophobicity.ts <- ggplot(data = as.data.frame(dat$pathwayProfileTimeSeries),
                                   aes(x = t,
                                       y = s,
                                       fill = pfHydrophobicity)) +
  geom_tile() +
  scale_fill_distiller(palette = "BrBG", 
                       direction = -1,
                       name = expression(paste(H~~bgroup("(", a.u.,")"))),
                       limits = c(-max(abs(dat$pathwayProfileTimeSeries$pfHydrophobicity)),
                                  max(abs(dat$pathwayProfileTimeSeries$pfHydrophobicity))),
                       guide = guide_colourbar(barheight = 10,
                                               barwidth = 1.3)) +
  scale_x_continuous(expand = c(0, 0),
                     name = expression(paste(t~~bgroup("(", ps,")")))) +
  scale_y_continuous(expand = c(0, 0),
                     name = expression(paste(s~~bgroup("(", nm,")")))) +
  ggtitle("Hydrophobicity Profile over Time") +
  theme_heatmap

ggsave("time_series_hydrophobicity_profile.png", 
       plt.pf.hydrophobicity.ts,
       width = plot.width,
       height = plot.height,
       dpi = opt$dpi)
