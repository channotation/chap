#!/usr/bin/env Rscript

# CHAP - The Channel Annotation Package
#
# Copyright (c) 2016 - 2019 Gianni Klesse, Shanlin Rao, Mark S. P. Sansom, and
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


# CONFIGURATION ================================================================

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

# get command line options:
option_list = list(
  make_option(
      c("--input", "-i"),
      action = "store",
      default = "output.json",
      type = "character",
      help = "CHAP output file in JSON format."
  ),
  make_option(
     c("--output", "-o"),
     action = "store",
     default = "profiles.csv",
     type = "character",
     help = "CHAP profiles file in CSV format."
  )
)
opt = parse_args(OptionParser(option_list=option_list))


# FORMAT CONVERSION ============================================================

# load chap output data:
message(
    paste0("loading CHAP data from ", opt$i, " ... "), appendLF = FALSE
)
chap_data <- fromJSON(readLines(opt$i, n = 1), flatten = FALSE)
message("done")

# write to CSV file:
message(
    paste0("writing pathway profiles to  ", opt$o, " ... "), appendLF = FALSE
)
write.csv2(
    chap_data$pathwayProfile,
    opt$o,
    row.names = FALSE
)
message("done")
