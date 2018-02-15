---
title: Molecular Graphics in PyMOL
permalink: /docs/molecular_graphics_pymol/
---

[WOBJ-spec]: http://www.fileformat.info/format/wavefrontobj/egff.htm
[PyMOL-main]: https://pymol.org/2/
[PyMOL-cgo]: https://pymolwiki.org/index.php/Category:CGO


CHAP writes a triangle mesh representation of the computed pore surface to a [Wavefront OBJ file][WOBJ-spec]. Unfortunately, [PyMOL][PyMOL-main] does not natively support the import of Wavefront OBJ meshes, which is why CHAP comes bundled with Python code that supports reading OBJ files and rendering them as [compiled graphics objects][PyMOL-cgo]. CHAP also provides a turnkey Python script for visualising the permeation pathway together with the channel protein.


## Loading and Displaying OBJ Files

The code for loading and displaying OBJ files can be found in the file `wobj.py`, which is located under `chap/scripts/visualisation/PyMOL/` in your CHAP install directory. Assuming that you copied `wobj.py` to your working directory, it can be imported into PyMOL via

```python
import wobj.py as wobj
```

and exposes the functions `import_wobj(filename)` and `draw_wobj(obj, groupname)`. 

The first of these can be used to load an OBJ file, which will be returned as a dictionary of compiled graphics objects (CGO), where each entry represents one group of faces loaded from the OBJ file, with the group names used as keys. The second function can subsequently be used to actually draw the pathway. So to display the permeation pathway coloured by solvent density, you would write:

```python
# load mesh data from file:
obj = wobj.import_wobj("output.obj")

# display surface:
wobj.draw_wobj(obj, groupname = "avg_density")
```

To see which groupnames are available, simply type `obj.keys()`. Note that the groupname argument to `draw_obj("output.obj")` is optional - if no groupname is specified, all groups will be drawn simultaneously and you can use the PyMOL GUI to hide and show the different surfaces.


## Turnkey Visualisation Script 

Under `chap/scripts/visualisation/PyMOL/` you can also find the Python script `visualise_pathway.py`, which will load both the `output.pdb` and `output.obj` files produced by CHAP to create a cartoon representation of the channel with the pore-lining and pore-facing residues additionally shown as sticks. This script will then use the functions defined in `wobj.py` to add the permeation pathway surface to the molecular visualisation. Therefore, if you want to use this script, you need to make sure that `wobj.py` is located either in your working directory or in the same location as `visualise_pathway.py`.

Assuming that you copied both scripts to your working directory, you can run the script either from within PyMOL by typing

```python
run visualise_pathway.py
```

or by passing it to the PyMOL start-up command:

```bash
pymol -r visualise_pathway.py --
```

Note that the two dashes are necessary for PyMOL to properly pass arguments to `visualise_pathway.py` even if no arguments are explicitly given.

By default, `visualise_pathway.py` assumes that the input files are named `output.pdb` and output `output.obj` and will produce one surface for each group in the OBJ file. When invoking the script from the command line, you can change this behaviour by passing named arguments after the double dash:

```bash
pymol -r visualise_pathway.py -- -structure custom.pdb -surface custom.obj -property avg_density
```
