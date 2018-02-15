---
title: Hydrophobicity Scales
permalink: /docs/hydrophobicity_scales/
---

[JSON-spec]: https://www.json.org/


CHAP associates each amino acid residue in a channel structure with a hydrophobicity value and uses this information to calculate a continuous hydrophobicity profile for the channel's permeation pathway as well as to produce hydrophobicity-coloured pore visualisations. In order to do this, it relies on published amino acid hydrophobicity scales. This page is intended to document the available hydrophobicity scales and to describe how custom scales can be added by a user.


## Hydrophobicity Scales in CHAP

The hydrophobicity scale data used by CHAP can be found in the `chap/share/data/hydrophobicity/` directory, with each scale contained in an individual [JSON][JSON-spec] file. Each file contains not only the scale data itself, but also a reference to the publication from which the data was taken. For convenience, these references are repeated in the table below: 

Scale Name 		| Reference  | DOI or Weblink
--- | --- | ---
`hessa_2005`            | Hessa, Tara, et al. "Recognition of transmembrane helices by the endoplasmic reticulum translocon." Nature 433.7024 (2005): 377-381.	| [10.1038/nature03216](https://doi.org/10.1038/nature03216)
`kyte_doolittle_1982`	| Kyte, Jack, and Russell F. Doolittle. "A simple method for displaying the hydropathic character of a protein." Journal of molecular biology 157.1 (1982): 105-132. |	[10.1016/0022-2836(82)90515-0](https://doi.org/10.1016/0022-2836(82)90515-0)
`memprotmd`		|  Newport, Thomas D. et al. "MemProtMD: A Database of Membrane Proteins Embedded in Lipid Bilayers" |  [www.memprotmd.bioch.oc.ac.uk/stats](http://memprotmd.bioch.ox.ac.uk/stats)
`monera_1995`		| Monera, Oscar D., et al. "Relationship of sidechain hydrophobicity and α‐helical propensity on the stability of the single‐stranded amphipathic α‐helix." Journal of peptide science 1.5 (1995): 319-329.	| [10.1002/psc.310010507](https://doi.org/10.1002/psc.310010507)
`moon_2015`		| Moon, C. Preston, and Karen G. Fleming. "Side-chain hydrophobicity scale derived from transmembrane protein folding into lipid bilayers." Proceedings of the National Academy of Sciences 108.25 (2011): 10174-10177. |	[10.1073/pnas.1103979108](https://doi.org/10.1073/pnas.1103979108)
`wimley_white_1996`	| Wimley, William C., and Stephen H. White. "Experimentally determined hydrophobicity scale for proteins at membrane interfaces." Nature Structural & Molecular Biology 3.10 (1996): 842-848.	| [10.1038/nsb1096-842](https://doi.org/10.1038/nsb1096-842) 
`zhu_2016`		|  Zhu, Chongqin, et al. "Characterizing hydrophobicity of amino acid side chains in a protein environment via measuring contact angle of a water nanodroplet on planar peptide network." Proceedings of the National Academy of Sciences (2016): 201616138. |	[10.1073/pnas.1616138113](https://doi.org/10.1073/pnas.1616138113)

In order to make the different scales comparable, hydrophobicity values have been linearly scaled to the interval [-1,1], with +1 being the most hydrophobic, -1 the most hydrophilic, and 0 a neutral residue. Where hydrophobicities have been scaled, the JSON file will also contain the originally published values in the `original_hydrophobicity` record.

Users can control which hydrophobicity scale is applied by setting the `-hydrophob-database` flag. By default, the `wimley_white_1996` scale is used. Note that all available scales only provide hydrophobicity data for amino acid residues. If the pathway-forming group (specified with the `-sel-pathway` flag) contains any other residues, you either have to provide a fallback hydrophobicity value with the `-hydrophob-fallback` flag or provide a custom hydrophobicity scale as described in the following section.


## Adding Custom Hydrophobicity Scales

In order to provide a custom hydrophobicity scale to CHAP, the `-hydrophob-database` flag needs to be set to `user` and the `-hydrophob-json` flag needs to point to a JSON file containing the scale data. 

The easiest way to create a hydrophobicity scale JSON file will be to adapt one of the existing files in `chap/share/data/hydrophobicity/`. Note that the only required record in this file is an object named `hydrophobicity` that contains an array of residue name and hydrophobicity value pairs as in the example below:

```json
{
    "hydrophobicity": [
        {"resname": "ALA", "hydrophobicity": -0.13692946},
        {"resname": "ARG", "hydrophobicity": -0.41493776},
        {"resname": "ASN", "hydrophobicity": -0.17842324},
        {"resname": "ASP", "hydrophobicity": -1.00000000},
        {"resname": "CYS", "hydrophobicity": -0.09128631},
        {"resname": "GLN", "hydrophobicity": -0.07883817},
        {"resname": "GLU", "hydrophobicity": -0.66804979},
        {"resname": "GLY", "hydrophobicity": -0.47302905},
        {"resname": "HIS", "hydrophobicity": -0.56846473},
        {"resname": "ILE", "hydrophobicity":  0.33609959},
        {"resname": "LEU", "hydrophobicity":  0.33609959},
        {"resname": "LYS", "hydrophobicity": -0.75103734},
        {"resname": "MET", "hydrophobicity":  0.18257261},
        {"resname": "PHE", "hydrophobicity":  0.24066390},
        {"resname": "PRO", "hydrophobicity":  0.12863071},
        {"resname": "SER", "hydrophobicity": -0.13692946},
        {"resname": "THR", "hydrophobicity": -0.04564315},
        {"resname": "TRP", "hydrophobicity":  0.09958506},
        {"resname": "TYR", "hydrophobicity": -0.09543568},
        {"resname": "VAL", "hydrophobicity":  0.21991701}
    ]
}
``` 

For consistency, you should rescale your hydrophobicity values to the interval [-1, 1] as described in the foregoing section.

