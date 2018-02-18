---
title: Annotation Example
permalink: /docs/annotation_example/
---

CHAP can be used with both PDB structures and molecular dynamics trajectories. The directory `examples` directory of the CHAP download contains example data for both cases and the following two sections will give a brief overview of how to run CHAP.


## Running CHAP on Structures

Under `examples/example-01/` you can find the file `4pirtm.pdb`, which contains the transmembrane domain of a serotonin receptor (PDB: 4PIR) channel embedded in a POPC bilayer and solvated in a sodium chloride solution. To compute a radius profile for this channel structure, type

```
chap -f 4pirtm.pdb -s 4pirtm.pdb
```

where the PDB file is used as both input structure (`-f`) and input topology (`-s`). You will be prompted to select the group of atoms forming the permeation pathway. Only atoms in this group will be considered in the probe-based pathway finding step of CHAP. Usually you will want to select the protein group, so press `1` and press enter. Note that you can also specify the pathway selection directly from the command line by typing:

```
chap -f 4pirtm.pdb -s 4pirtm.pdb -sel-pathway 1
```

Once you pressed enter, CHAP will load data from `4pirtm.pdb`, compute the permeation pathway, identify pore-lining residues, and calculate a pathway hydrophobicity profile accordingly. In order to visualise the results, run

```
./path/to/chap/scripts/plotting/Python/chap_plot_pathway_profile.py
```

which will generate four plots, one for pore radius, pore hydrophobicity, solvent density, and solvent free energy each (a corresponding R script is also available). In addition to the profiles themselves, the radius and hydrophobicity plots will also show a scatter plot of the pore-lining residues.

As you may notice, the density and energy profiles are simply flat lines. The next section will illustrate how to obtain meaningful results for these profiles.


## Running CHAP on Trajectories

In the directory `examples/example-02` you can find the files `4pirtm.xtc` and `4pirtm.tpr`. These are a trajectory and topology files for the same serotonin receptor membrane system as discussed in the previous section. The trajectory only contains 10 frames covering 10 ns from a much longer simulation, but this will be sufficient to illustrate how CHAP works on trajectories. Run CHAP on this data by typing:

```
chap -f 4pirtm.xtc -s 4pirtm.tpr -sel-pathway 1 -sel-solvent 16
```

Note that in addition to the pathway selection the above command also specifies the solvent selection as all water particles in the system. CHAP will never prompt you for a solvent selection and if you do not specify it explicitly, CHAP will only compute radius and hydrophobicity profiles, but not solvent density and free energy profiles.

When running the above command, CHAP will go through the trajectory frame by frame and run its pathway-finding and solvent-density-estimation algorithms on each frame. Subsequently, CHAP will calculate profile time-averages and further summary statistics (such as the local standard deviation). Again, you can run

```
./path/to/chap/scripts/plotting/Python/chap_plot_pathway_profile.py
```

to generate pathway profile plots. In this case, the solvent density and free energy profiles should no longer be flat lines and will exhibit a clear dip and peak respectively. You may also notice a grey shaded area around each profile. The dark gray shade corresponds to the one standard deviation range of the profile and the light grey shaded area corresponds to the minimum-maximum range for that profile over all trajectory frames. The dashed blue line in the solvent density profile plot corresponds to the literature value for bulk water number density of ca 33.37 nm<sup>-3</sup>.

You can also add a visualisation of the pathway to a molecular graphic of the protein. To do this, copy the files from `path/to/chap/scripts/visualisation/VMD` to the working directory and then run VMD by typing:

```
vmd -e visualise_pathway.tcl -args output.pdb output.obj avg_energy
```

This will generate a representation of the protein as ribbons with pore-lining residues represented by their van der Waals spheres and a surface representation of the permeation pathway colour-coded by the solvent free energy.


## Next Steps ##

The above is a very brief illustration of how CHAP works and what sort of output it generates. You may want to take a look at the [Visualising Results](http://www.channotation.org/docs/plotting_python/) chapter for a more detailed description of how to plot CHAP results in both R and Python and for how to generate molecular graphics in both VMD and PyMOL. You may also be interested in reading the chapter on [Understanding Output](http://www.channotation.org/docs/overview_output/) for a comprehensive overview of the data made available by CHAP and the various output files. The chapter on [Running CHAP](http://www.channotation.org/docs/chap_workflow/) provides an in-depth discussion of how CHAP works and what parameters can be used to tweak its performance. Finally, the chapter on [Data Sources](http://www.channotation.org/docs/vdw_radii/) lists the academic publications from which the van der Waals radius and residue hydrophobicity data used in CHAP is drawn.

