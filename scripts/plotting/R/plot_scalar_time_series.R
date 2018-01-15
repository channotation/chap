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
dat <- fromJSON(readLines(filename, n = 1), flatten = FALSE)


################################################################################
# TIME SERIES PLOTS
################################################################################

# individual plot size:
plot.width <- unit(plot.width.cm, "cm")
plot.height <- unit(plot.height.cm, "cm")

# smoother length:
t.smooth <- 5000 # time scale of smoother
span.smooth <- t.smooth/abs(diff(range(as.data.frame(dat$pathwayScalarTimeSeries)$t)))


# Length
#-------------------------------------------------------------------------------

plt.length <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                     aes(x = t, 
                         y = length)) +
  geom_line(colour = brewer.pal(9, "Greys")[6]) +
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

ggsave("time_series_pathway_length_over_time.png", 
       plt.length,
       width = plot.width,
       height = plot.height,
       dpi = 1200)


  # Volume
#-------------------------------------------------------------------------------

plt.volume <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                     aes(x = t, 
                         y = volume)) +
  geom_line(colour = brewer.pal(9, "Greys")[6]) +
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
       dpi = 1200)


# Number of Solvent Particles in Pathway
#-------------------------------------------------------------------------------

plt.solvent.number <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                             aes(x = t, 
                                 y = numPathway)) +
  geom_line(colour = brewer.pal(9, "Greys")[6]) +
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

ggsave("time_series_avg_solvent_number_density.png", 
       plt.solvent.number,
       width = plot.width,
       height = plot.height,
       dpi = 1200)


# Average Solvent Number Density in Pathway
#-------------------------------------------------------------------------------

plt.avg.solvent.number.density <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                                         aes(x = t, 
                                             y = numPathway/volume)) +
  geom_line(colour = brewer.pal(9, "Greys")[6]) +
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
       dpi = 1200)


# Minimal Pore Radius
#-------------------------------------------------------------------------------

plt.min.pore.radius <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                              aes(x = t, 
                                  y = minRadius)) +
  geom_line(colour = brewer.pal(9, "Greys")[6]) +
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
       dpi = 1200)


# Minimal Solvent Number Density
#-------------------------------------------------------------------------------

plt.min.solvent.density <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                                  aes(x = t, 
                                      y = minSolventDensity)) +
  geom_line(colour = brewer.pal(9, "Greys")[6]) +
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
       dpi = 1200)


# Location of Minimal Pore Radius
#-------------------------------------------------------------------------------

plt.argmin.pore.radius <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                                 aes(x = t, 
                                     y = argMinRadius)) +
  geom_point(colour = brewer.pal(9, "Greys")[6]) +
  scale_x_continuous(expand = c(0,0)) + 
  xlab(expression(paste(t~~bgroup("(", ps,")")))) +
  ylab(expression(paste(arg*min(R, s)~~bgroup("(", nm,")")))) +
  ggtitle("Location of Minimal Pore Radius over Time") +
  theme_chap

ggsave("time_series_argmin_pore_radius.png", 
       plt.argmin.pore.radius,
       width = plot.width,
       height = plot.height,
       dpi = 1200)


# Location of Minimal Solvent Number Density
#-------------------------------------------------------------------------------

plt.argmin.solvent.density <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                                     aes(x = t, 
                                         y = argMinSolventDensity)) +
  geom_point(colour = brewer.pal(9, "Greys")[6]) +
  scale_x_continuous(expand = c(0,0)) + 
  xlab(expression(paste(t~~bgroup("(", ps,")")))) +
  ylab(expression(paste(arg*min(n, s)~~bgroup("(", nm^{-3},")")))) +
  ggtitle("Location of Minimal Density over Time") +
  theme_chap

ggsave("time_series_argmin_solvent_number_density.png", 
       plt.argmin.solvent.density,
       width = plot.width,
       height = plot.height,
       dpi = 1200)


# Combined Plot
#-------------------------------------------------------------------------------

plt.timeseries <- grid.arrange(plt.length, 
                               plt.volume,
                               plt.solvent.number,
                               plt.avg.solvent.number.density,
                               plt.min.pore.radius,
                               plt.min.solvent.density,
                               plt.argmin.pore.radius,
                               plt.argmin.solvent.density,
                               ncol = 2)

ggsave("time_series_scalar.png", 
       plt.timeseries, 
       width = unit(2*plot.width.cm, "cm"), 
       height = unit(4*plot.height.cm, "cm"),
       dpi = 300)


