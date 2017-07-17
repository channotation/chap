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
cd chap
mkdir build
cd build
~~~


## Usage ##


```
#!bash

./chap -chan-dir-vec 0 0 1 -init-probe-pos 0 0 0 -max-free-dist 1 -pf-method optim-direction -probe-radius 1 -probe-step 0.2 -sa-conv-tol 1e-3 -sa-cooling-fac 0.9 -sa-cost-samples 10 -sa-init-temp 100 -sa-max-cool 1000 -sa-random-seed 15011991 -sa-step 0.1 -f pr.gro -n pr.ndx
```

## Visualisation ##

CHAP will write PDB files containing the positions and radii of the probe spheres, where the radii are contained in the "beta" column. To visualise path:

1. Load PDB file into VMD.
2. In Tcl console issue `set sel [atomselect top "name PORE"]` to select all pore pseudo-particles.
3. Then issue `$sel set radius [$sel get beta]` to get correct radii.
