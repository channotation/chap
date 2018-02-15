---
title: Colour Palettes
permalink: /docs/colour_palettes/
---


[colorbrewer]: http://colorbrewer2.org/
[ggplot2]: http://ggplot2.org/
[Matplotlib]: https://matplotlib.org/

For colour-coding pathway surfaces that are written to Wavefront OBJ and MTL files, CHAP uses a set of colour palettes, which can be found under `chap/share/data/palettes/default.json`. This file specifies one colour scale for each pathway property. These colour scales have been created by lab-space interpolation of the discrete palettes from the [colorbrewer][colorbrewer] project, which are natively available in many common plotting packages such as R's [ggplot2][ggplot2] or Python's [Matplotlib][Matplotlib]. In particular, the following scales are used:

Property 		| Palette	| Remark
--- 			| --- 		| ---
`radius`		| `YlOrBr`	| Darker hues indicate smaller radii. 
`avg_radius`		| `YlOrBr`	| Darker hues indicate smaller radii.
`avg_pl_hydrophobicity`	| `BrBG`	| Hydrophobic is brown, hydrophilic is blue-green, neutral is white.
`avg_pf_hydrophobicity` | `BrBG`	| Hydrophobic is brown, hydrophilic is blue-green, neutral is white.
`avg_density`		| `Blues`	| Darker hues indicate larger density.
`avg_energy`		| `PuRd`	| Darker hues indicate higher energy.

The palettes used are chosen to be colourblind-safe, intuitive, and widely available, but please keep in mind that their appearance in molecular graphics will also depend on factors over which CHAP has no influence, such as the lighting of the scene or the material properties.

