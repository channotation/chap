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
plot.width.cm <- 5.8
plot.height.cm <- plot.width.cm * 0.75

# plot appearance:
theme_heatmap <- theme(text = element_text(size = 13),
                       axis.text.x = element_text(size = 12,
                                                  colour = "black"),
                       axis.text.y = element_text(size = 12,
                                                  colour = "black"),
                       panel.border = element_rect(colour = "black", 
                                                   fill=NA, 
                                                   size=1)) 
  
  
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
       dpi = 1200)


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
       dpi = 1200)


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
       dpi = 1200)


# Combined Plot
#-------------------------------------------------------------------------------

plt.profile.ts <- grid.arrange(plt.radius.profile.ts, 
                               plt.density.ts,
                               plt.pf.hydrophobicity.ts,
                               ncol = 1)

ggsave("time_series_pathway_profiles.png", 
       plt.profile.ts, 
       width = unit(1*plot.width.cm, "cm"), 
       height = unit(3*plot.height.cm, "cm"),
       dpi = 400)
