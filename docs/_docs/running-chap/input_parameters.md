---
title: Input Parameters
permalink: /docs/input_parameters/
---

[HOLE]: http://www.holeprogram.org/


The following is a comprehensive list of input parameters to CHAP. Note that you can get a similar overview including available options and default values in CHAP's command line help by typing `chap -h`.


## Input Files

Controls what data is loaded into CHAP. Note that a topology is always required, so if you want to run CHAP on an individual structure use `chap -f structure.pdb -s structure.pdb` and for a trajectory use `chap -f trajectory.xtc -s topology.tpr`.

`-f`    |   Input trajectory or single configuration.
`-s`    |   Input topology.
`-n`    |   Index file for custom index groups.


## Time Options

These control the time window from which frames are loaded and the spacing of windows. Note that if you want to estimate energies you should discount the beginning of your trajectory as equilibration time.

`-b`    |   First frame (in picoseconds) to read from trajectory.
`-e`    |   Last frame (in picoseconds) to read from trajectory.
`-dt`   |   Only use frame if t MOD dt == first time (in picoseconds).
`-tu`   |   Unit for time values: fs, ps, ns, us, ms, s.


## Unused Gromacs Options

These are added by `libgromacs` per default, but are unused in CHAP.

`-fgroup`   | unused
`-xvg`      | unused
`-sf`       | unused
`-selrpos`  | unused
`-seltype`  | unused


## Pathway and Solvent Selections

These two options tell CHAP what constitutes the permeation pathway and solvent.

Only atoms in the selection specified by `-sel-pathway` will be considered in the path-finding algorithm, all other atoms will be ignored. Usually, this flag will be set to the `Protein` group. This option is mandatory and CHAP will ask for it on start-up if it has not been specified on the command line. 

Atoms which are part of the selection specified by `-sel-solvent` will be considered in the density estimation step. Usually, this flag will be set to the `Water` group. This flag is optional and if no solvent selection is specified, the density profile in the output data will simply be zero. 

`-sel-pathway`  | Reference group that defines the permeation pathway.
`-sel-solvent`  | Group of small particles to calculate density of.


## Output Options

Control the name of the output files and allow for some tweaking of the pathway surface written to the output OBJ and MTL files.

---                 | ---
`-out-filename`     |   File name for output files without file extension. 
`-out-num-points`   |   Number of spatial sample points that are written to the JSON output file.
`-out-extrap-dist`  |   Extrapolation distance beyond the pathway endpoints for both JSON and OBJ output.
`-out-grid-dist`    |   Controls the sampling distance of vertices on the pathway surface which are subsequently interpolated to yield a smooth surface. Very small values may yield visual artefacts.
`-out-vis-tweak`    |    Visual tweaking factor that controls the smoothness of the pathway surface in the OBJ output. Varies between -1 and 1 (exclusively), where larger values result in a smoother surface. Negative values may result in visualisation artefacts.
`-[no]out-detailed` |   If true, CHAP will write detailed per-frame information to a newline-delimited JSON file including original probe positions and spline parameters. This is mostly useful for debugging.


## Pathway-Finding Options

These parameters control how CHAP determines the permeation pathway through the group of atoms specified by `-sel-pathway`.

The default pathway-finding method `inplane_optim` uses an algorithm similar to that used in the [HOLE][HOLE] programme, where a spherical probe is squeezed through the van der Waals spheres of the pathway-forming atoms. The specific van der Waals radii used in the procedure can be set using the `-pf-vdwr-database`, `-pf-vdwr-fallback`, and `-pf-vdwr-json` flags.

The initial probe position is normally chosen to be the centre of mass (COM) of the pathway-forming atoms, but this behaviour can be changed through the `-pf-sel-ipp` (COM of a different group of atoms is used) and `-pf-init-probe-pos` (initial probe position is specified explicitly) flags. The probe is then moved by `-pf-probe-step` in the direction of the channel direction vector specified with `-pf-chan-dir-vec`. Note that by default this vector points in the Cartesian z-direction and for most ion channels it is more sensible to align the channel protein appropriately than to adjust the channel direction vector.

The probe motion is stopped if either a pathway radius larger than `-pf-max-free-dist` is encountered or the probe has already moved by `-pf-max-probe-steps` steps. The point at which this happens will be considered the pathway endpoint and the probe is then moved in the opposite direction of `-pf-chan-dir-vec` to find the other pathway endpoint.

Alternatively, the `-pf-method` flag can be set to `cylindrical` if the above method fails to find the correct pathway. In this case, the permeation pathway will be a cylindrical volume centred around the initial probe position and extending `-pf-max-probe-steps` times `-pf-probe-step` in either direction along the axis specified by `-pf-chan-dir-vec`. Note that in general the `cylindrical` method will not produce an accurate radius profile for the permeation pathway and consequently the solvent density profile will not take into account a variation of free space along the pathway.

`-pf-method`            |   Pathway-finding method.
`-pf-vdwr-database`     |   Database of van der Waals radii to be used in pathway finding.
`-pf-vdwr-fallback`     |   Fallback van der Waals radius for atoms that are not listed in van der Waals radius database.
`-pf-vdwr-json`         |   JSON file with user-defined van der Waals radii. Will be ignored unless `-pf-vdwr-database` is set to `user`.
`-pf-align-method`      |   Method for aligning pathway coordinates across time steps.
`-pf-probe-step`        |   Step length for probe movement.
`-pf-max-free-dist`     |   Maximum radius of pore. The point at which this radius is reached marks the endpoint of the pathway.
`-pf-max-probe-steps`   |   Maximum number of steps the probe is moved in either direction.
`-pf-sel-ipp`           |   Selection of atoms whose COM will be used as initial probe position. If not set, the selection specified with `-sel-pathway` will be used.
`-pf-init-probe-pos`    |   Initial position of probe in probe-based pore finding algorithms. If set explicitly, it will overwrite the COM-based initial position set with `-sel-ipp`.
`-pf-chan-dir-vec`      |   Channel direction vector. Will be normalised to unit vector internally.
`-pf-cutoff`            |   Cutoff distance for spatial searches in pathway-finding algorithm. A value of zero or less means no cutoff is applied. If unset, a cutoff is determined automatically.


## Optimisation Parameters used in Pathway Finding

The probe-based pathway finding-algorithm outlined above uses two subsequent optimisation procedures for finding the position which maximises the probe radius: A global optimisation procedure based on simulated annealing and a local optimisation procedure based on the Nelder-Mead simplex method. The parameters below can be used to tweak these optimisation methods. Note simulated annealing is turned off by default.

`-sa-seed`          |   Seed used in pseudo random number generation for simulated annealing. If not set explicitly, a random seed is used.
`-sa-max-iter`      |   Number of cooling iterations in one simulated annealing run.
`-sa-init-temp`     |   Simulated annealing initial temperature.
`-sa-cooling-fac`   |   Simulated annealing cooling factor.
`-sa-step`          |   Step length factor used in candidate generation.
`-nm-max-iter`      |   Number of Nelder-Mead simplex iterations.
`-nm-init-shift`    |   Distance of vertices in initial Nelder-Mead simplex.


## Pathway-Mapping Parameters

In order to determine which residues are lining the permeation pathway, CHAP calculates the distance between each residue's centre of geometry and the pathway centre line. If this distance is smaller than the local pathway radius plus `-pm-pl-margin` the residue is considered to be *pore-lining*. If in addition the residue's Cα atom is further from the centre line than its centre of geometry, a residue will also be considered *pore-facing*, i.e. its side chain is pointing towards the pore.

Note that for determining whether a residue is pore-facing, CHAP requires each residue in the pathway forming group to contain a Cα atom, which will usually be the case for protein channels. However, if the selection specified with `-sel-pathway` contains non-amino-acid residues (e.g. carbon nanotubes, DNA nanopores, or protein channels where lipids form part of the pathway), this condition will not be met. In this case, CHAP will assume that no residue is pore-facing and only properties calculated from pore-lining residues will be meaningful.

In order to allow custom reference positions for the determination of pore-facing residues, the `-pm-pf-sel` flag can be used. The centre of geometry of the subset of each residue specified by `-pm-pf-sel` will then take the place of the Cα position and a residue will be considered pore-lining, if this centre of geometry is further away from the centre line of the pathway than the centre of geometry of the whole residue. 

`-pm-pl-margin`     |   Margin for determining pathway-lining residues.
`-pm-pf-sel`        |   Selection string determining the centre of geometry group for assessing if a residue is pore-facing. 


## Density Estimation Parameters

In order to determine the solvent density along the permeation pathway, CHAP first maps the COM position of all residues in the `-sel-solvent` selection onto the pathway centre line. Subsequently, it uses the method specified with the `-de-method` flag to estimate the one-dimensional probability density of residue positions.

By default, a kernel density estimator with an automatically determined bandwidth is used, but the bandwidth can also be set explicitly with the `-de-bandwidth` flag or fine-tuned with the `-de-bw-scale` flag. If a histogram is used for density estimation, the `-de-res` flag can be used to specify the histogram bin width; for a kernel estimator this parameter determines the spacing of evaluation points.

`-de-method`        |   Method used for estimating the probability density of the solvent particles along the permeation pathway.
`-de-res`           |   Spatial resolution of the density estimator. In case of a histogram, this is the bin width, in case of a kernel density estimator, this is the spacing of the evaluation points.
`-de-bandwidth`     |   Bandwidth for the kernel density estimator. Ignored for other methods. If negative or zero, bandwidth will be determined automatically.
`-de-bw-scale`      |   Scaling factor for the band width. Useful to set a bandwidth relative to the automatically determined value.
`-de-eval-cutoff`   |   Evaluation range cutoff for kernel density estimator in multiples of bandwidth. Ignored for other methods. Ensures that the density falls off smoothly to zero outside the data range.


## Hydrophobicity Parameters

In addition to radius and solvent density profiles, CHAP also computes a hydrophobicity profile. This is accomplished by kernel smoothing of hydrophobicity values associated with the pore-lining residues. The hydrophobicity associated with each residue can be controlled through the `-hydrophob-database`, `-hydrophob-fallback`, and `-hydrophob-json` flags. The amount of smoothing can be controlled with `-hydrophob-bandwidth`, where larger values will generate a smoother profile.

`-hydrophob-database`   |   Database of hydrophobicity scale for pore-forming residues.
`-hydrophob-fallback`   |   Fallback hydrophobicity for residues in the pathway-defining group.
`-hydrophob-json`       |   JSON file with user-defined hydrophobicity scale. Will be ignored unless `-hydrophob-database` is set to `user`.
`-hydrophob-bandwidth`  |   Bandwidth for hydrophobicity kernel.

