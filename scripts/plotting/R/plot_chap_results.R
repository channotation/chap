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
  








