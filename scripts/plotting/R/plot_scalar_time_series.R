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
  geom_line() +
  geom_smooth(method = "loess",
              se = FALSE,
              colour = "gold",
              span = span.smooth,
              n = 2048,
              size = 1.2) +
  xlab(expression(paste(t~~bgroup("(", ps,")")))) +
  ylab(expression(paste(L~~bgroup("(", nm,")")))) +
  ggtitle("Permeation Pathway Length over Time") +
  theme_chap

ggsave("pathway_length_over_time.pdf", 
       plt.length,
       width = plot.width,
       height = plot.height)


  # Volume
#-------------------------------------------------------------------------------

plt.volume <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                     aes(x = t, 
                         y = volume)) +
  geom_line() +
  geom_smooth(method = "loess",
              se = FALSE,
              colour = "gold",
              span = span.smooth,
              n = 2048,
              size = 1.2) +
  xlab(expression(paste(t~~bgroup("(", ps,")")))) +
  ylab(expression(paste(V~~bgroup("(", nm^{3},")")))) +
  ggtitle("Permeation Pathway Volume over Time") +
  theme_chap

ggsave("pathway_volume_over_time.pdf", 
       plt.volume,
       width = plot.width,
       height = plot.height)


# Number of Solvent Particles in Pathway
#-------------------------------------------------------------------------------

plt.solvent.number <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                             aes(x = t, 
                                 y = numPathway)) +
  geom_line() +
  geom_smooth(method = "loess",
              se = FALSE,
              colour = "gold",
              span = span.smooth,
              n = 2048,
              size = 1.2) +
  xlab(expression(paste(t~~bgroup("(", ps,")")))) +
  ylab(expression(paste(N))) +
  ggtitle("Number of Solvent Particles in Pathway over Time") +
  theme_chap

ggsave("avg_solvent_number_density_over_time.pdf", 
       plt.solvent.number,
       width = plot.width,
       height = plot.height)


# Average Solvent Number Density in Pathway
#-------------------------------------------------------------------------------

plt.avg.solvent.number.density <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                                         aes(x = t, 
                                             y = numPathway/volume)) +
  geom_line() +
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
  xlab(expression(paste(t~~bgroup("(", ps,")")))) +
  ylab(expression(paste(N/V~~bgroup("(", nm^{-3},")")))) +
  ggtitle("Average Solvent Number Density in Pathway over Time") +
  theme_chap

ggsave("avg_solvent_number_density_over_time.pdf", 
       plt.avg.solvent.number.density,
       width = plot.width,
       height = plot.height)


# Minimal Pore Radius
#-------------------------------------------------------------------------------

plt.min.pore.radius <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                              aes(x = t, 
                                  y = minRadius)) +
  geom_line() +
  geom_smooth(method = "loess",
              se = FALSE,
              colour = "gold",
              span = span.smooth,
              n = 2048,
              size = 1.2) +
  xlab(expression(paste(t~~bgroup("(", ps,")")))) +
  ylab(expression(paste(min(R, s)~~bgroup("(", nm,")")))) +
  ggtitle("Minimal Pore Radius over Time") +
  theme_chap

ggsave("min_pore_radius_over_time.pdf", 
       plt.min.pore.radius,
       width = plot.width,
       height = plot.height)


# Minimal Solvent Number Density
#-------------------------------------------------------------------------------

plt.min.solvent.density <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                                  aes(x = t, 
                                      y = minSolventDensity)) +
  geom_line() +
  geom_smooth(method = "loess",
              se = FALSE,
              colour = "gold",
              span = span.smooth,
              n = 2048,
              size = 1.2) +
  xlab(expression(paste(t~~bgroup("(", ps,")")))) +
  ylab(expression(paste(min(n, s)~~bgroup("(", nm^{-3},")")))) +
  ggtitle("Minimal Solvent Number Density over Time") +
  theme_chap

ggsave("min_solvent_number_density_over_time.pdf", 
       plt.min.solvent.density,
       width = plot.width,
       height = plot.height)


# Location of Minimal Pore Radius
#-------------------------------------------------------------------------------

plt.argmin.pore.radius <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                                 aes(x = t, 
                                     y = argMinRadius)) +
  geom_point() +
  xlab(expression(paste(t~~bgroup("(", ps,")")))) +
  ylab(expression(paste(arg*min(R, s)~~bgroup("(", nm,")")))) +
  ggtitle("Location of Minimal Pore Radius over Time") +
  theme_chap

ggsave("argmin_pore_radius_over_time.pdf", 
       plt.argmin.pore.radius,
       width = plot.width,
       height = plot.height)


# Location of Minimal Solvent Number Density
#-------------------------------------------------------------------------------

plt.argmin.solvent.density <- ggplot(data = as.data.frame(dat$pathwayScalarTimeSeries),
                                     aes(x = t, 
                                         y = argMinSolventDensity)) +
  geom_point() +
  xlab(expression(paste(t~~bgroup("(", ps,")")))) +
  ylab(expression(paste(arg*min(n, s)~~bgroup("(", nm^{-3},")")))) +
  ggtitle("Location of Minimal Solvent Number Density over Time") +
  theme_chap

ggsave("argmin_solvent_number_density_over_time.pdf", 
       plt.argmin.solvent.density,
       width = plot.width,
       height = plot.height)


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

ggsave("time_series.pdf", 
       plt.timeseries, 
       width = unit(2*plot.width.cm, "cm"), 
       height = unit(4*plot.height.cm, "cm"))


