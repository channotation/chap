# README #

## Installation ##

1. From `build` directory, issue `cmake ..`.
2. Then issue `make -j 12` from same directory.
3. From `bin` directory issue `./runAllTests` to make sure compilation went well.


## Usage ##


```
#!bash

./chap -chan-dir-vec 0 0 1 -init-probe-pos 0 0 0 -max-free-dist 1 -pf-method hole -probe-radius 1 -probe-step 0.2 -sa-conv-tol 1e-3 -sa-cooling-fac 0.9 -sa-cost-samples 10 -sa-init-temp 100 -sa-max-cool 1000 -sa-random-seed 15011991 -sa-step 0.1 -f pr.gro -n pr.ndx
```