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

plt.radius.profile <- ggplot(data = as.data.frame(dat$residueSummary)[dat$residueSummary$poreFacing$mean >= 0.1 &
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
  xlab(expression(paste(s~~bgroup("(",nm,")")))) +
  ylab(expression(paste(R~~bgroup("(",nm,")")))) +
  ggtitle("Time-Averaged Pathway Radius Profile") +
  theme_chap

ggsave("time_averaged_radius_profile.pdf", 
       plt.radius.profile,
       width = plot.width,
       height = plot.height)


# Radius Profile with Residue Positions
#-------------------------------------------------------------------------------

plt.pf.residues <- ggplot(data = as.data.frame(dat$residueSummary)[dat$residueSummary$poreFacing$mean >= 0.1 &
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
                 colour = name),
             size = 2) +
  xlab(expression(paste(s~~bgroup("(",nm,")")))) +
  ylab(expression(paste(R~~bgroup("(",nm,")")))) +
  ggtitle("Time-Averaged Pore-Facing Residue Positions") +
  theme_chap

ggsave("time_averaged_residue_positions.pdf", 
       plt.pf.residues,
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
                            plt.pf.residues,
                            plt.solvent.number.density,
                            plt.free.energy,
                            ncol = 2)

ggsave("pathway_profiles.pdf", 
       plt.profile, 
       width = unit(2*plot.width.cm, "cm"), 
       height = unit(2*plot.height.cm, "cm"))

