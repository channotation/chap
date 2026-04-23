#!/bin/bash

mkdir build
cd build
# Ensure installed binaries can find shared libs in conda-style prefixes without requiring LD_LIBRARY_PATH.
cmake -DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX} \
	-DCMAKE_INSTALL_RPATH="\$ORIGIN/../lib;\$ORIGIN/../lib.AVX2_256;\$ORIGIN/../../lib;\$ORIGIN/../../lib.AVX2_256" \
	-DCMAKE_BUILD_WITH_INSTALL_RPATH=OFF \
	-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=OFF \
	..
make
make check
make install
cp ${CONDA_PREFIX}/chap/bin/* ${CONDA_PREFIX}/bin
