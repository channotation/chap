# README #

## Installation ##

1. From `build` directory, issue `cmake ..`.
2. Then issue `make -j 12` from same directory.
3. Issue `make install`.
4. From `bin` directory issue `./runAllTests` to make sure compilation went well.


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