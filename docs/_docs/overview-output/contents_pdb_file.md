---
title: Contents of the PDB Output File
permalink: /docs/contents_pdb_file/
---

The `output.pdb` written by CHAP reproduces the structure contained in the
topology (the input passed to CHAP with the `-s` flag), but overwrites the
occupancy and temperature factor fields of each `ATOM` record with values
indicating whether an atom belongs to a residue that is considered pore-lining
and pore-facing respectively.

To be precise, both fields contain a number between zero and one that indicates
the fraction of time a residue was considered pore-lining or pore-facing
respectively, where a value of one means that a residue has been pore-lining (or
pore-facing) for every frame in the input trajectory, while a value of zero
indicates that a residue has never been found to be pore-lining (or
pore-facing).  As the pore-lining and pore-facing attributes are calculated on a
per-residues basis, this number will be identical for all atoms of the same
residue.

