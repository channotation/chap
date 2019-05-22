#!/usr/bin/env Rscript

# CHAP - The Channel Annotation Package
# 
# Copyright (c) 2016 - 201 Gianni Klesse, Shanlin Rao, Mark S. P. Sansom, and 
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


# Configuration ================================================================

if (!require(optparse)) {
  install.packages("optparse", repos = "http://cran.us.r-project.org")
  library(optparse)
}
if (!require(jsonlite)) {
  install.packages("jsonlite", repos = "http://cran.us.r-project.org")
  library(jsonlite)
}
if (!require(ggplot2)) {
  install.packages("ggplot2", repos = "http://cran.us.r-project.org")
  library(ggplot2)
}
if (!require(scales)) {
  install.packages("scales", repos = "http://cran.us.r-project.org")
  library(scales)
}
if (!require(gridExtra)) {
  install.packages("gridExtra", repos = "http://cran.us.r-project.org")
  library(gridExtra)
}
if (!require(extrafont)) {
  install.packages("extrafont", repos = "http://cran.us.r-project.org")
  library(extrafont)
  font_import()
}

option_list <- list(
  make_option(c("--filename"),
              action = "store",
              default = "output.json",
              type = "character",
              help = "Name of the input file."),
  make_option(c("--energies"),
              action = "store",
              default = "heuristic_grid.json",
              type = "character",
              help = "Name of the SVM model file.")
)
opt <- parse_args(OptionParser(option_list = option_list))


# Data read-in =================================================================

chap <- fromJSON(readLines(opt$filename, n = 1), flatten = FALSE)
grid <- fromJSON(readLines(opt$energies, n = 1), flatten = FALSE)


# Plot appearance ==============================================================

theme_heuristic <- theme_bw() +
  theme(panel.border = element_rect(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 13, family = "Ubuntu"),
        axis.title = element_text(size = 14, family = "Ubuntu"),
        plot.title = element_text(hjust = 0.5, size = 15, family = "Ubuntu"))


# Radius and hydrophobicity profiles ===========================================

radius.profile <- ggplot(data = as.data.frame(chap$pathwayProfile),
                         aes(x = s, y = radiusMean)) +
  geom_hline(yintercept = 0.14,
             colour = "#006080",
             size = 0.6,
             linetype = "dotted") +
  geom_line() +
  scale_x_continuous(breaks = seq(-30, 30, 2), expand = c(0.002, 0)) +
  scale_y_continuous(breaks = seq(0, 2, 0.2), expand = c(0.002, 0)) +
  coord_cartesian(ylim = c(0, 1.1)) +
  xlab("s (nm)") + ylab("Radius (nm)") +
  ggtitle("") +
  theme_heuristic +
  theme(plot.margin = unit(c(0.2, 0.2, 0, 0.2), "cm"))

hydrophobicity.profile <- ggplot(as.data.frame(chap$pathwayProfile),
                                 aes(x = s, y = pfHydrophobicityMean)) +
  geom_line() +
  scale_x_continuous(breaks = seq(-30, 30, 2), expand = c(0.002, 0)) +
  scale_y_continuous(expand = c(0.002, 0)) +
  coord_cartesian(ylim = c(-1, 0.4)) +
  xlab("s (nm)") + ylab("Hydrophobicity") +
  ggtitle("") +
  theme_heuristic +
  theme(plot.margin = unit(c(0, 0.2, 0.2, 0.2), "cm"))


# Prepare data frame for pore-lining residues ==================================

narrow <- chap$residueSummary$poreRadius$mean < 0.7
pore_facing <- chap$residueSummary$poreFacing$mean > 0.5

df <- c()

for (residue in chap$residueSummary$id[narrow & pore_facing]) {

  nearest <- which.min(abs(chap$pathwayProfile$s
                           - chap$residueSummary$s$mean[chap$residueSummary$id == residue]))

  temp <- data.frame(name = chap$residueSummary$name[chap$residueSummary$id == residue],
                     hydrophobicity = chap$pathwayProfile$pfHydrophobicityMean[nearest],
                     radius = chap$pathwayProfile$radiusMean[nearest])

  df <- rbind(df, temp)
}


# Heuristic prediction for structure ===========================================

for (i in 1:nrow(df)) {
  df$pred[i] <- grid$energy[which.min( (grid$hydrophobicity - df$hydrophobicity[i])^2
                                      + (grid$radius - df$radius[i])^2)]
}

mapped <- ggplot() +
  geom_contour(data = grid,
               aes(x = hydrophobicity,
                   y = radius,
                   z = energy,
                   colour = ..level..),
               breaks = c(2.6),
               linetype = "dashed",
               size = 0.5) +
  geom_point(data = df,
             aes(x = hydrophobicity,
                 y = radius),
             fill = ifelse(df$pred > 2.6, "#9c4842", "#ffffff"),
             shape = 22,
             size = 3,
             stroke = 0.6,
             colour = "grey20") +
  scale_colour_gradient(low = "#466a83",
                        high = "#ffbc98",
                        limits = c(0, 10.4),
                        oob = squish,
                        guide = FALSE) +
  scale_y_continuous(expand = c(0.002, 0)) +
  scale_x_continuous(expand = c(0.002, 0)) +
  coord_cartesian(ylim = c(0.09, 0.7)) +
  theme_heuristic +
  theme(plot.margin = unit(c(0.2, 0.7, 0.2, 0.2), "cm")) +
  ggtitle("Prediction result") +
  labs(y = "Pore radius (nm)", x = "Hydrophobicity")

divide <- ggplot_build(mapped)[[1]][[1]]
matrix <- matrix(c(divide$x, divide$y), ncol = 2)
points <- ggplot_build(mapped)[[1]][[2]][ggplot_build(mapped)[[1]][[2]]$fill == "#9c4842", ]
distances <- points[, c(2, 1)]

if (nrow(points) != 0) {
  for (i in 1:nrow(points)) {
    distances$shortest[i] <-
      min(sqrt(rowSums(sweep(matrix, 2, as.numeric(points[i, c(2, 1)]), "-")^2)))
  }
}

SSD <- round(sum(distances$shortest), digits = 2)

labelled <- mapped +
  annotate("text",
           x = 0,
           y = 0.64,
           label = paste("(n = ", nrow(distances), ")", sep = ""),
           family = "Ubuntu",
           size = 4,
           colour = "#33526f") +
  annotate("text",
           x = 0,
           y = 0.58,
           label = "heuristic score:",
           family = "Ubuntu",
           size = 5,
           colour = "#33526f") +
  annotate("text",
           x = 0,
           y = 0.53,
           label = SSD,
           family = "Ubuntu",
           size = 5,
           colour = "#33526f")


# Arrange and save plots =======================================================

profiles <- grid.arrange(labelled,
                         arrangeGrob(radius.profile,
                                     hydrophobicity.profile,
                                     nrow = 2),
                         ncol = 2,
                         widths = c(7, 5))

ggsave("heuristic_prediction.png",
       profiles,
       dpi = 600,
       height = 3.7,
       width = 6.4)
