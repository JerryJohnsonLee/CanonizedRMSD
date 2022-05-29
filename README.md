# A Bifunctional Canonization for Efficient Minimal Root-mean-sqaure Deviation (RMSD) Calculation and CIP Stereochemistry Identification 

This repository provides a tool that can be used to canonize molecules and provide CIP stereochemistry tags for atoms in the molecule.  The branching tiebreaking process implemented in the algorithm provides the minimal RMSD between two molecules.

## Usage notes

usage: `CanonizedRMSD.py [-h] [-s] [-m] [-i] [-q] [-r] file1 file2`

to calculate the RMSD of two molecules after canonizing them.

supported file types: `.mol | .sdf | .rxn | .mol2 | .ml2 | .pdb`

positional arguments: `file1 file2`

optional arguments:

<pre>
-h, --help show this help message and exit

-s, --save save intermediate results

-m, --mapping output atom mapping relationship with two molecules

-i, --ignore_isomerism ignore geometric and stereometric isomerism when canonizing

-a, --no_alignment do not apply molecule alignment by Kabsch algorithm or QCP algorithm when calculating RMSD

-q, --use QCP algorithm instead of Kabsch algorithm

-r, --remove H atoms
</pre>