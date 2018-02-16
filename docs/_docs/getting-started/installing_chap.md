---
title: Installing CHAP
permalink: /docs/installing_chap/
---

[CMake]: https://cmake.org/
[Boost]: http://www.boost.org/
[CBLAS]: http://www.netlib.org/blas/
[LAPACKE]: http://www.netlib.org/lapack/lapacke.html
[Gromacs]: http://www.gromacs.org/
[Gromacs-install]: http://manual.gromacs.org/documentation/ 
[GCC]: https://gcc.gnu.org/
[GTest]: https://github.com/google/googletest


CHAP is written in C++ and is developed and tested on Ubuntu Linux. It is expected to be portable to other Unix-like operating systems, although this has not been tested. The following installation guide will assume that you are using Ubuntu and that you are reasonably experienced in using the command line.


## Prerequisites ##

Prior to installing CHAP, make sure that you have the following libraries and tools installed:

1. The [CMake][CMake] tool in version 3.2 or higher. This will typically be available through your system's package manager. For example, on Ubuntu you can install CMake by typing `sudo apt-get install cmake`. CMake is used to check the availability of libraries and compilers on your system and will ensure that CHAP is installed properly.
2. A C++ compiler that supports the `C++11` standard. A popular choice is the [GNU Compiler Collection][GCC], which on Ubuntu can be obtained by typing `sudo apt-get install gcc`
3. The [Boost][Boost] C++ libraries, which on Ubuntu can be installed using `sudo apt-get install libboost-all-dev`. Boost algorithms are used in CHAP to solve some root finding and optimisation problems.
4. The CBLAS and LAPACKE linear algebra libraries. On Ubuntu, the easiest way to obtain these is by typing `sudo apt-get install libblas-dev liblapacke-dev libatlas-base-dev`. The linear algebra libraries are used in CHAP's spline interpolation.
5. The `libgromacs` library of the [Gromacs][Gromacs] molecular dynamics engine in version 2016 or higher. Comprehensive installation instructions for Gromacs can be found [here][Gromacs-install].
Please note that for using Gromacs as a library, the underlying FFTW libray 
may **not** be installed automatically, i.e. you need to set
`-DGMX_BUILD_OWN_FFTW=OFF` when running CMake during the Gromacs 
installation.

CHAP also depends on [RapidJSON](http://rapidjson.org/), but this is included as a header-only library, and on [GTest][GTest], but this is downloaded and installed automatically by CMake, so you don't need to do anything about either of these (you will however need Internet access when installing CHAP).


## Compiling and Installing CHAP  ##

To install CHAP, unpack the source, create a `build` directory parallel to the source tree and from there run `cmake`, `make`, `make check`, and `make install`.

```bash
cd chap
mkdir build
cd build
cmake ..
make
make check
sudo make install
```

CMake will automatically find all dependencies (and inform you of missing ones), `make` will compile the code (you can use the `make -j` flag to speed this up on multicore machines), `make check` runs a suite of unit tests, and `make install` will place the binary in `/usr/local/chap/bin` (so you will need sudo rights for this last step). To check if CHAP has been installed properly, you can type

```bash
chap -h
```

which should bring up an online help for using CHAP.

