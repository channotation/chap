################################################################################
# SETTINGS (CHANGE FILENAME HERE)
################################################################################

# libraries:
library(ggplot2)   # plotting
library(ggrepel)   # improved labels
library(gridExtra) # plot panel arrangement
library(jsonlite)  # parsing JSON files

# name of data file:
setwd("/sansom/s117/scro2967/year-2/misc/faradary-discussion-paper/figures/pore_profiles")
filename <- "output.json"

# plot output parameters:
plot.width.cm <- 21.0/2
plot.height.cm <- 29.7/4

# plot appearance:
theme_chap <- theme_grey(base_size = 18)
theme_heatmap <- theme(text = element_text(size = 13),
                       axis.text.x = element_text(size = 12,
                                                  colour = "black",
                                                  margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
                       axis.text.y = element_text(size = 12,
                                                  colour = "black"),
                       panel.border = element_rect(colour = "black", fill=NA, size=1)) 
  
  
nm.to.ang <- 10


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


# Radius Profile Time Series
#-------------------------------------------------------------------------------

ggplot(data = as.data.frame(dat$pathwayProfileTimeSeries),
       aes(x = t,
           y = s,
           fill = radius)) +
  geom_tile() +
  scale_fill_distiller(palette = "YlOrRd",
                       direction = -1,
                       name = expression(paste(R~~bgroup("(", nm,")")))) +
  scale_x_continuous(expand = c(0, 0),
                     name = expression(paste(t~~bgroup("(", ps,")")))) +
  scale_y_continuous(expand = c(0, 0),
                     name = expression(paste(s~~bgroup("(", nm,")")))) +
  theme_heatmap


# Solvent Density Profile Time Series
#-------------------------------------------------------------------------------

ggplot(data = as.data.frame(dat$pathwayProfileTimeSeries),
       aes(x = t,
           y = s,
           fill = density)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues",
                       direction = 1,
                       name = expression(paste(n~~bgroup("(", nm^{-3},")")))) +
  scale_x_continuous(expand = c(0, 0),
                     name = expression(paste(t~~bgroup("(", ps,")")))) +
  scale_y_continuous(expand = c(0, 0),
                     name = expression(paste(s~~bgroup("(", nm,")")))) +
  theme_heatmap
  
  
# Pore-lining Hydrophobicity Profile Time Series
#-------------------------------------------------------------------------------

ggplot(data = as.data.frame(dat$pathwayProfileTimeSeries),
       aes(x = t,
           y = s,
           fill = plHydrophobicity)) +
  geom_tile() +
  scale_fill_distiller(palette = "BrBG", 
                       direction = -1,
                       name = expression(paste(H~~bgroup("(", a.u.,")"))),
                       limits = c(-max(abs(dat$pathwayProfileTimeSeries$plHydrophobicity)),
                                  max(abs(dat$pathwayProfileTimeSeries$plHydrophobicity)))) +
  scale_x_continuous(expand = c(0, 0),
                     name = expression(paste(t~~bgroup("(", ps,")")))) +
  scale_y_continuous(expand = c(0, 0),
                     name = expression(paste(s~~bgroup("(", nm,")")))) +
  theme_heatmap



# Pore-facing Profile Time Series
#-------------------------------------------------------------------------------

ggplot(data = as.data.frame(dat$pathwayProfileTimeSeries),
       aes(x = t,
           y = s,
           fill = pfHydrophobicity)) +
  geom_tile() +
  scale_fill_distiller(palette = "BrBG", 
                       direction = -1,
                       name = expression(paste(H~~bgroup("(", a.u.,")"))),
                       limits = c(-max(abs(dat$pathwayProfileTimeSeries$pfHydrophobicity)),
                                  max(abs(dat$pathwayProfileTimeSeries$pfHydrophobicity)))) +
  scale_x_continuous(expand = c(0, 0),
                     name = expression(paste(t~~bgroup("(", ps,")")))) +
  scale_y_continuous(expand = c(0, 0),
                     name = expression(paste(s~~bgroup("(", nm,")")))) +
  theme_heatmap
