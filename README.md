# README #

## Installation ##

### Prerequisites ###

Prior to installing CHAP, make sure that you have the following libraries and tools installed:

1. The [CMake](https://cmake.org/) tool in version 2.8.8 or higher. This will 
typically available through your system's package manager. For example, on 
Ubuntu you can install CMake by typing `sudp apt-get install cmake`.
2. The BLAS and LAPACK linear algebra libraries. This will typically already be
installed on your system and are otherwise available through your package 
manager. On Ubuntu, type `sudo apt-get install libblas-dev liblapack-dev` to
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
`make1 will compile the code (you can use the `make -j` flag to speed this up 
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

For a complete list of all options and a brief online helpyou can type 
`chap -h`.


## Visualising Results ##

### Output File Structure ###

The results file returned by CHAP is organised as a newline delimted JSON file
which means that each line is a valid JSON object that can be read and parsed
individually. The first line contains results aggregated over all trajectory
frames and each subsequent line contains information for one trajectory frame
(in the appropriate order). As an end user you will typically only need to look
at the first line and under `chap/scripts/plotting/` you can find ready made
scripts for visualising the results.

## Plotting CHAP Results in R ##

An example script for plotting CHAP results in R can be found under 
`chap/scripts/plotting/R/plot_chap_results.R`. 


## Plotting CHAP Results in Python ##

Python scripts for plotting CHAP results are currently under development. 
Please be patient.

