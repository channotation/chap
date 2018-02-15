---
title: Molecular Graphics in VMD
permalink: /docs/molecular_graphics_vmd/
---

[WOBJ-spec]: http://www.fileformat.info/format/wavefrontobj/egff.htm
[VMD-main]: http://www.ks.uiuc.edu/Research/vmd/
[VMD-draw]: http://www.ks.uiuc.edu/Research/vmd/current/ug/node127.html


CHAP writes a triangle mesh representation of the computed pore surface to a [Wavefront OBJ file][WOBJ-spec]. Unfortunately, [VMD][VMD-main] does not natively support the import of Wavefront OBJ meshes, which is why CHAP comes bundled with Tcl code that supports reading and visualising OBJ files using [VMD's draw command][VMD-draw]. CHAP also provides a turnkey Tcl script for visualising the permeation pathway together with the channel protein.


## Loading and Displaying OBJ Files

The code for loading and displaying OBJ data can by found in the file `wobj.tcl`, which is located under `chap/scripts/visualisation/VMD/` in your CHAP install directory. Assuming that you copied `wobj.tcl` to your working directory, you can load it into VMD by opening the Tk Console (under `Extensions > Tk Console` in the VMD menu) and typing:

```tcl
source wobj.tcl
```

This will make the functions `import_wobj` and `draw_wobj` available. The first of these can be used to load data from an OBJ file and the second function can be used to draw either the entire mesh or individual groups of faces. Both functions are contained within a namespace called `WOBJ`, so to load data from `output.obj`, you need to type:

```tcl
set obj [WOBJ::import_wobj "output.obj"]
```

You can now display the permeation pathway coloured by e.g. solvent density by typing:

```
WOBJ::draw_wobj $obj "avg_density"
```

Note that if the second argument to `draw_wobj` is a string that does not match any of the group names in `output.obj` the pore will be drawn in solid colour.

As VMD's `draw` command adds elements to the molecule as a whole, it is unfortunately not possible to add different surface colours as individual representations that can be switched off or on through the GUI and pathway visualisation has to be controlled manually through VMD's Tk console.


## Turnkey Visualisation Script 

To allow quick visualisation of CHAP results in VMD, the directory `chap/scripts/visualisation/VMD/` also contains the Tcl script `visualise_pathway.tcl`, which loads data from both `output.pdb` and `output.obj` to visualise the channel protein together with permeation pathway. Internally, this script uses functions defined in `wobj.tcl`, so you need to make sure that `wobj.tcl` is located either in your working directory or in the same directory as `visualise_pathway.tcl`.

The `visualise_pathway.tcl` script can be started either from within VMD by typing

```tcl
source visualise_pathway.tcl
```

in the Tk console or directly from the command line as part of the VMD invocation by typing:

```bash
vmd -e visualise_pathway.tcl
```

If the second approach is used, `wobj.tcl` needs to be in the current working directory.

By default, `visualise_pathway.tcl` assumes that the data files are called `output.obj` and `output.pdb`. In order to use custom file names, you need to set the variables `FILE_PORE_SURFACE` and `FILE_STRUCTURE` before calling `source visualise_pathway.tcl`. Similarly, if you want the pathway surface coloured by a specific property, you need to set the variable `PROPERTY` accordingly. For example, if your data files are named `abcd.pdb` and `pathway.obj` and you want the pathway surfaced coloured by solvent density, you need to type the following:

```tcl
set FILE_STRUCTURE "abcd.pdb"
set FILE_PORE_SURFACE "pathway.obj"
set PROPERTY "avg_density"
source visualise_pathway.tcl
```

You can also pass arguments to `visualise_pathway.tcl` if you are loading it at VMD start-up, by using the `-args` flag to set the PDB file name, OBJ file name, and property name (in this order):

```tcl
vmd -e visualise_pathway.tcl -args abcd.pdb pathway.obj avg_density
```

Note that in this case all three arguments need to be set explicitly.
