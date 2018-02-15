---
title: Contents of the JSON Output File
permalink: /docs/contents_json_file/
---

[JSON-spec]: https://www.json.org/
[JQ]: https://stedolan.github.io/jq/


The `output.json` file is the main output file of CHAP and contains all data
except for what is needed in 3D molecular visualisation. It abides by the
specification of the [JSON file format][JSON-spec] and can easily be read into
any common programming language for visualisation or further processing. 

By default, CHAP writes its output to a *minified* JSON file, which means that
when you open `output.json` in a text editor, you will find all data in one line
and without any white spaces. This is done to save disk space and speed up the 
reading and writing process, but does not result in a human-readable file. You
can use a tool such as [jq][JQ] to *prettify* the JSON file. If you run

```sh
jq . output.json >> prettified.json
```

the resulting `prettified.json` file will contain nicely reformatted JSON with
just one record per line. This is helpful to get an overview of how
`output.json` is structured, but for most use cases it will be more appropriate
to read the minified JSON file into the scripting language of your choice and to
process the data therein. 

On the highest level, `output.json` contains six JSON objects, which are
summarised in the table below:

Object Name                  | Summary
---                          | ---
`reproducibilityInformation` | Information on the CHAP version used and the parameters with which it was called.
`pathwaySummary`             | Summary statistics for scalar-valued properties of the channel, such as its volume or minimal radius.
`pathwayProfile`             | Profiles of properties varying along the channel, such as local radius or hydrophobicity.
`pathwayScalarTimeSeries`    | Time series for scalar-valued channel properties.
`pathwayProfileTimeSeries`   | Time series for properties varying along the channel.
`residueSummary`             | Summary statistics on various residue properties.

The remainder of this chapter will provide a detailed description of the
information contained in each of these JSON objects.


## Reproducibility Information

The reproducibility information is intended to help you retrace the steps
undertaken to obtain a result. It contains information on the version of CHAP 
that generated a given data file and the command line call used to invoke chap.
In practice, it may look like this:

```json
{
  "reproducibilityInformation": {
    "version": {
      "string": "0.7.0",
      "major": "0",
      "minor": "7",
      "patch": "0",
      "gitHash": "e8766d6074820abbe93981cdf486f8a2a8327975"
    },
    "commandLine": "chap -f pr.xtc -s pr.tpr -sel-pathway 1 -sel-solvent 16"
  }
}
```

Here `version` is an object that contains the version number, which is
generally written as a dot-separated string of the major, minor, and patch
version number. The `gitHash` entry is an automatically generated serial number
that changes every time the underlying code of CHAP changes. The `commandLine`
entry is simply a string that can be copied to re-execute CHAP with the same 
parameters as the run that generated `output.json`.


## Pathway Summary

The pathway summary provides information on the scalar-valued variables that 
describe the permeation pathway as a whole. Each of these quantities is given
as a set of summary statistics (in particular the minimum, maximum, mean, 
standard deviation, and variance) that describe how the variable fluctuates over 
time. Inside `output.json` the pathway summary may look like this: 

```json
{
  "pathwaySummary": {
    "argMinRadius": {
      "min": 7.898504257202148,
      "max": 8.483589172363281,
      "mean": 8.265358924865723,
      "sd": 0.15247607231140137,
      "var": 0.023248951882123947
    },
    "minRadius": { 
		...
    },
    "length": {
		...
    },
    "volume": {
		...
    },
    "numPathway": {
		...
    },
    "numSample": {
		...
    },
    "argMinSolventDensity": {
		...
    },
    "minSolventDensity": {
		...
    },
    "bandWidth": {
		...
    }
  }
}
```

For brevity, the summary statistics have been omitted for all but the first 
variable in the above example. The following table gives a brief description of 
each of the variables contained within `pathwaySummary`:

Variable | Description
--- | ---
`length` 				| The length of the permeation pathway measured as the distance between the two openings of the pathway pore (i.e. the points where the radius reaches `-pf-max-free-dist`). 
`volume` 				| The volume of the permeation pathway measured as the integral over cross-sectional area between the two openings of the pathway pore.
`minRadius`				| The minimum radius of the pore between the two openings of the pathway pore.
`argMinRadius`			| The location of the minimum radius along the pathway centre line.
`numPathway`			| The number of solvent particles inside the permeation pathway.
`numSample`				| The number of solvent particles inside the sample used to calculate solvent number density.
`minSolventDensity`		| The minimum solvent number density between the two openings of the pore.
`argMinSolventDensity`	| The location of the minimum solvent number density along the pathway centre line.
`bandWidth`				| The bandwidth used in the kernel density estimate of the solvent probability density.


## Pathway Profile

The `pathwayProfile` object contains time-averaged data on variables that vary
along the length of the permeation pathway. It is laid out as an object of
multiple same-length JSON arrays (you can think of this as a data table). The 
first JSON array contains the pathway coordinate `s`; all other arrays contain
summary statistics of pathway properties evaluated at the given value of `s`. 
In particular, the minimum, maximum, mean, and standard deviation of each 
variable are available and the array name is a composition of the variable name
and the summary statistic:

```json
{
  "pathwayProfile":{
    "s": [...],
    "radiusMin": [...],
    "radiusMax": [...],
    "radiusMean": [...],
    "radiusSd": [...],
    "densityMin": [...],
    "densityMax": [...],
    "densityMean": [...],
    "densitySd": [...],
    "energyMin": [...],
    "energyMax": [...],
    "energyMean": [...],
    "energySd": [...],
    "...": [...]
  }
}
```
Note that for brevity not all properties are explicitly listed in the above 
example and that for further clarity the number contained within each array 
have been omitted. The following table gives a comprehensive overview of all
pathway properties in `pathwayProfile`. The summary statistic suffix is omitted
here.

Variable 			| Description
--- 				| ---
`s`					| The arc length along the centre line of the permeation pathway. By default, the zero is at the pathway forming selection's centre of mass, but this may change based on the (`-pf-align-method`, `-pf-sel-ipp`, and `-pf-init-probe-pos` flags).
`radius`			| Radius of the permeation pathway. Defined as the maximum radius a spherical probe located on the centre line can have without overlapping any pore atom's van der Waals radius. Will be primarily influenced by the `-pf-*` flags.
`plHydrophobicity`  | Hydrophobicity of the permeation pathway due to pore-lining residues. This is calculated via kernel smoothing of the hydrophobicities associated with the pore-lining residues and is influenced by the `-hydrophob-*` and `-pm-*` flags.
`pfHydrophobicity`	| Hydrophobicity of the permeation pathway due to pore-facing residues. This is calculated via kernel smoothing of the hydrophobicities associated with the pore-facing residues and is influenced by the `-hydrophob-*` and `-pm-*` flags. 
`density`			| Number density of solvent particles along the pathway. If no solvent particle selection is given, this array just contains an arbitrary constant. Primarily influenced by `-de-*` flags.
`energy`			| Free energy profile of solvent particles. Only meaningful for data generated from sufficently long and well-equilibrated trajectories. If no solvent particle selection is given, this array just contains an arbitrary constant. Calculated directly from number density and therefore influenced by the same flags.

Note that the range and granularity of `s` values for which profile data is 
written to `output.json` can be controlled with the `-out-extrap-dist` and
`-out-num-points` flags.


## Pathway Scalar Time Series

The `pathwayScalarTimeSeries` object contains time series data arrays for all 
scalar-valued pathway properties in `pathwaySummary`. In addition, it contains
an array of time stamps and is overall of the following form:

```json
{
  "pathwayScalarTimeSeries": {
    "t": [...],
    "argMinRadius": [...],
    "minRadius": [...],
    "length": [...],
    "volume": [...],
    "numPathway": [...],
    "numSample": [...],
    "argminSolventDensity": [...],
    "minSolventDensity": [...],
    "bandWidth": [...]
  }
}
```


## Pathway Profile Time Series

The `pathwayProfileTimeSeries` object contains time series for all properties 
for which time-averaged data exists in the `pathwayProfile` object with the 
exception of energy, which is only meaningfully defined as time-averaged
property. It contains one JSON array for each pathway property as well as two 
additional JSON arrays containing the time stamps and spatial coordinates:

```
{
  "pathwayProfileTimeSeries": {
    "t": [...],
    "s": [...],
    "radius": [...],
    "density": [...],
    "plHydrophobicity": [...],
    "pfHydrophobicity": [...]
  }
}
```

The `pathwayProfileTimeSeries` object can therefore be thought of as a 
long-format data table, with time being the leading dimension (i.e. the value
of `t` is repeated as many times as there are different values of `s`).


## Residue Summary

The `residueInformation` object contains data on a per-residue basis and in 
particular the following variables:

Variable 			| Description
--- 				| ---
`id`				| Residue ID as in input topology.
`name`				| Residue name as in input topology.
`chain`				| Residue chain as in input topology.
`hydrophobicity`	| Residue hydrophobicity as in hydrophobicity database (see `-hydrophob-database` flag).
`s`					| Summary statistics for residue COM position along the pathway centre line.
`rho`				| Summary statistics for residue COM distance from centre line.
`phi`				| Not currently meaningfully defined.
`poreLining`		| Summary statistics for pore-lining attribute.
`poreFacing`		| Summary statistics for pore-facing attribute.
`poreRadius`		| Summary statistics for pore radius at residue position.
`solventDensity`	| Summary statistics for solvent density at residue position.
`x`					| Summary statistics for residue COM Cartesian x-coordinate.
`y`					| Summary statistics for residue COM Cartesian y-coordinate.
`z`					| Summary statistics for residue COM Cartesian z-coordinate.

For all but the first four variables, which are essentially non-varying database
values, the minimum, maximum, mean, standard deviation, and variance (over time)
are available.


## Units and Further Notes

By default CHAP output contains the following units:

Dimension 		| Unit
--- 			| ---
time			| ps (can be changed with `-tu` flag)
length			| nm (but Angstrom used in OBJ files)
number density  | nm<sup>-3</sup>
energy			| k<sub>B</sub>T (where T is given in Kelvin)

Please also keep in mind the following additional notes:

 * CHAP defines the variance and standard deviation to be *zero* in cases where fewer than two samples are present. While a more mathematically intuitive approach would define these quantities as infinite in these cases, setting both quantities to zero allows reusing plot scripts with error bars.
 * Due to a technical limitation of the JSON format, CHAP output can not contain infinities (which may occur as energy values where the density drops to zero). Wherever infinities do occur these are written to output as the largest representable floating point number, with the sign being the same as the sign of the infinity.

