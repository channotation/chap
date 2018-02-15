---
title: Overview of Output Files
permalink: /docs/overview_output/
---

[WOBJ-spec]: http://www.fileformat.info/format/wavefrontobj/egff.htm
[WMTL-spec]: http://www.fileformat.info/format/material/
[JSON-spec]: https://www.json.org/
[JQ]: https://stedolan.github.io/jq/


CHAP will generally produce four different output files each time it is run. By
default, these are called `output.json`, `output.pdb`, `output.obj`, and
`output.mtl`, but the base name (i.e. the part before the respective file
extension) can be changed using the `-out-filename` flag. The JSON file contains
all information you need to plot profiles of the permeation pathway (e.g. its
radius along the centre line) and their evolution over time in a machine-readable 
format. The other three files contain the information that is needed to
visualise the pathway in a molecular visualisation system such as PyMOL or VMD.

The next few sections document the contents of the output files and are intended
as a comprehensive reference. If you are mainly interested in quickly
visualising your CHAP results, take a look at the menu on the left, where you
may find practical examples of how to plot CHAP results in Python and R and how
to create molecular visualisations of the permeation pathway in PyMOL and VMD.

