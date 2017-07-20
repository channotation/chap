# README #

## Installation ##

### Prerequisites ###

Prior to installing CHAP, make sure that you have the following libraries and tools installed:

1. The [CMake](https://cmake.org/) tool in version 2.8.8 or higher. This will 
typically available through your system's package manager. For example, on 
Ubuntu you can install CMake by typing `sudp apt-get install cmake`.
2. The BLAS and LAPACKE linear algebra libraries. This will typically already be
installed on your system and are otherwise available through your package 
manager. On Ubuntu, type `sudo apt-get install libblas-dev liblapacke-dev` to
install both.
3. The [Google Test](https://github.com/google/googletest) unit testing 
framework. This needs to be compiled from source using CMake.
4. The `libgromacs` library of the [Gromacs](http://www.gromacs.org/) molecular 
dynamics engine. Very detailed installation instructions for this can be found
[here](http://manual.gromacs.org/documentation/2016.3/install-guide/index.html).
Please note that for using Gromacs as a library, the underlying FFTW libray 
may **not** be installed automatically, i.e. you need to set
`-DGMX_BUILD_OWN_FFTW=OFF` to false when running CMake during the Gromacs 
installation.

CHAP also depends on [RapidJSON](http://rapidjson.org/), but this is included
as a header only library, so you don't need to do anything about this.


### Quick and Dirty Installation ###

To install CHAP, unpack the source, create a `build` directory parallel to the
source tree and from there run `cmake`, `make`, `make check`, and 
`make install`.

~~~
#!bash
cd chap
mkdir build
cd build
cmake ..
make
make check
sudo make install
~~~

CMake will automatically find all dependencies (and inform you of missing ones)
`make` will compile the code (you can use the `make -j` flag to speed this up 
on multicore machines to speed this up), `make check` runs a suite of unit 
tests, and `make install` will place the binary in `/usr/local/chap` (so you
need sudo rights for this last step.


## Usage ##

### Running CHAP ###

To run chap on a MD trajectory, you simply need to to provide the name of the
trajectory file and the corresponding topology (which is required to assign
a van-der-Waals radius to each particle when finding the permeation pathway):

```
#!bash
chap -f trajectory.xtc -s topology.tpr
```

This will create a file `output.json` which contains all results in newline
delimited [JSON](http://www.json.org/) format. The next section describes how
the results file is organised and how to visualise your results.

For a complete list of all options and a brief online help you can type 
`chap -h`.


## Visualising Results ##

### Output File Structure ###

The results file returned by CHAP is organised as a newline delimted JSON file
which means that each line is a valid JSON object that can be read and parsed
individually. The first line contains results aggregated over all trajectory
frames and each subsequent line contains information for one specific 
trajectory frame (in the appropriate order). As an end user you will typically 
only need to look at the first line and under `chap/scripts/plotting/` you can 
find ready made scripts for visualising the results.

### Plotting CHAP Results in R ###

Start by reading the first line of the CHAP output file into R and parse the 
results using a library such as 
[jsonlite](https://cran.r-project.org/web/packages/jsonlite/index.html)

~~~
library(jsonlite) # parsing JSON files
dat <- fromJSON(readLines(filename, n = 1), flatten = FALSE)
~~~

The now contains a named list that maintains the structure of the results file.

In `dat$pathSummary` you can find a summary of the overall properties of the
permeation pathway such as its length, volume, and minimum radius. All these
quantities aggregated over time and you are provided a set of summary 
statistics (minimum, maximum, mean, and standard deviation) that indicate how 
much these variables varied over time.

The object `dat$pathProfile` represents a table of pathway radius, solvent
density, and solvent free energy along the centre line of the pore (this object
can be coerced to a data frame). All three quantities are again aggregated 
over time and you are again provided with a set of summary statistics.

As an example, you can visualise the radius along the permation pathway and its
variation over time using the following script:
~~~
library(ggplot2)
ggplot(as.data.frame(dat$pathProfile),
       aes(x = s,
           y = radiusMean)) +
  geom_ribbon(alpha = 0.4,
              aes(ymin = radiusMean - radiusSd,
                  ymax = radiusMean + radiusSd)) +
  geom_ribbon(alpha = 0.2,
              aes(ymin = radiusMin,
                  ymax = radiusMax)) +
  geom_ribbon(alpha = 0.4) +
  geom_line()
~~~
This will generate a line graph of the time averaged radius with two ribbons
representing the one standard deviation interval and the range between minimum
and maximum radius. 

Lastly, in `dat$reproducibilityInformation` you can find some information on 
the CHAP version used to generate the data as well as the parameters that were
used. This is intended to help you to reproduce your results as the software
evolves.

An example script for plotting CHAP results in R can be found under 
`chap/scripts/plotting/R/plot_chap_results.R`. 


### Plotting CHAP Results in Python ###

Python scripts for plotting CHAP results are currently under development. 
Please be patient.


### Molecular Visualisation in VMD ###

Visualising CHAP pore surfaces in VMD is possible, but this component of the
program is still under development and contains some known and possibly 
more unknown bugs. Moreover, the workflow will likely change from what
is described below in the near future.

In addtion to the main results file in JSON format, CHAP produces two more
output files, `output.obj` and `res_mapping.dat`. These contain a pore 
surface mesh in 
[Wavefront OBJ](https://en.wikipedia.org/wiki/Wavefront_.obj_file) format
data identifying the pore lining and pore facing residues respectively.

Under `chap/scripts/visualisation/vmd` you can find the `pore_lining.tcl`
script, which allows you to visualise the pore together with your protein 
molecule in VMD. To do this you need to open the Tk Console in VMD and 
specify the names of three files:

~~~
set FILE_STRUCTURE structure.gro
set FILE_PORE_LINING res_mapping.dat
set FILE_PORE_SURFACE output.obj
~~~

You can then load the data by sourcing the `pore_lining.tcl` script:

~~~
source pore_lining.tcl
~~~

Note however that this script uses a itself sources the `wobj.tcl` script
(which contains an OBJ parser), which must be present in the same directory.


### Molecular Visualisation in PyMol ###

PyMol scripts for visualising CHAP results are currently under development.
Please be patient.
