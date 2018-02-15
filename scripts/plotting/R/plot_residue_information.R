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
# SETTINGS
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
# RESIDUE INFORMATION PLOTS
################################################################################


# Positions of Pore-Facing Residues
#-------------------------------------------------------------------------------

data <- as.data.frame(dat$residueSummary)
plt.residue.positions <- ggplot(data = data,
                                aes(x = s.mean,
                                    xmin = s.mean - s.sd,
                                    xmax = s.mean + s.sd,
                                    y = rho.mean,
                                    ymin = rho.mean - rho.sd,
                                    ymax = rho.mean + rho.sd,
                                    colour = poreFacing.mean)) +
  geom_point(size = 2.5) +
  scale_colour_distiller(palette = "PRGn",
                         name = expression(paste("pore facing"~~bgroup("(","%",")")))) +
  xlab(expression(paste(s~~bgroup("(",nm,")")))) +
  ylab(expression(paste(rho~~bgroup("(",nm,")")))) +
  ggtitle("Time-Averaged Residue Positions") +
  theme_chap 

ggsave("residue_positions.pdf", 
       plt.residue.positions,
       width = plot.width,
       height = plot.height)


# Hydrophobicity of Pore-Facing Residues
#-------------------------------------------------------------------------------

data <- as.data.frame(dat$residueSummary)
plt.hydrophobicity <- ggplot(data = data[data$poreFacing.mean > 0.1,],
                                aes(x = s.mean,
                                    xmin = s.mean - s.sd,
                                    xmax = s.mean + s.sd,
                                    y = rho.mean,
                                    ymin = rho.mean - rho.sd,
                                    ymax = rho.mean + rho.sd,
                                    colour = hydrophobicity)) +
  geom_point(size = 2.5) +
  geom_errorbar(size = 0.75) +
  geom_errorbarh(size = 0.75) +
  scale_colour_distiller(palette = "RdBu",
                         name = expression(paste("hydrophobicity  "))) +
  xlab(expression(paste(s~~bgroup("(",nm,")")))) +
  ylab(expression(paste(rho~~bgroup("(",nm,")")))) +
  ggtitle("Hydrophobicity of Pore-Facing-Residues") +
  theme_chap 
  
ggsave("residue_positions.pdf", 
       plt.residue.positions,
       width = plot.width,
       height = plot.height)


# Pore Radius at Pore-Facing Residues
#-------------------------------------------------------------------------------

data <- as.data.frame(dat$residueSummary)
plt.pore.radius.at.res <- ggplot(data = data[data$poreFacing.mean > 0.1,],
                                 aes(x = as.factor(id),
                                     y = poreRadius.mean, 
                                     ymin = poreRadius.mean - poreRadius.sd,
                                     ymax = poreRadius.mean + poreRadius.sd,
                                     colour = name)) +
  geom_point(size = 3) +
  geom_errorbar(size = 1.2) +  
  xlab("residue ID") +
  ylab(expression(paste(R~~bgroup("(",nm,")")))) +
  ggtitle("Time-Averaged Pore Radius at Pore-Facing Residues") +
  theme_chap +
  theme(axis.text.x = element_text(angle=45,
                                   size = 11))

ggsave("pathway_radius_at_pf_residues.pdf", 
       plt.pore.radius.at.res,
       width = plot.width,
       height = plot.height)


# Solvent Number Density at Pore-Facing Residues
#-------------------------------------------------------------------------------

data <- as.data.frame(dat$residueSummary)
plt.solvent.density.at.res <- ggplot(data = data[data$poreFacing.mean > 0.1,],
                                     aes(x = as.factor(id),
                                         y = solventDensity.mean, 
                                         ymin = solventDensity.mean - solventDensity.sd,
                                         ymax = solventDensity.mean + solventDensity.sd,
                                         colour = name)) +
  geom_point(size = 3) +
  geom_errorbar(size = 1.2) +  
  xlab("residue ID") +
  ylab(expression(paste(n~~bgroup("(",nm^{-3},")")))) +
  ggtitle("Time-Averaged Solvent Density at Pore-Facing Residues") +
  theme_chap +
  theme(axis.text.x = element_text(angle=45,
                                   size = 11))

ggsave("solvent_density_at_pf_residues.pdf", 
       plt.solvent.density.at.res,
       width = plot.width,
       height = plot.height)


# Combined Plot
#-------------------------------------------------------------------------------

plt.residues <- grid.arrange(plt.residue.positions,
                             plt.pore.radius.at.res,
                             plt.hydrophobicity,
                             plt.solvent.density.at.res,
                             ncol = 2)

ggsave("residue_summary.pdf", 
       plt.residues, 
       width = unit(2*plot.width.cm, "cm"), 
       height = unit(2*plot.height.cm, "cm"))
