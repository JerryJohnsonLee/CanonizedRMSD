# A Bifunctional Canonization for Efficient Minimal Root-mean-sqaure Deviation (RMSD) Calculation and CIP Stereochemistry Identification 

This repository provides a tool that can be used to canonize molecules and provide CIP stereochemistry tags for atoms in the molecule.  The branching tiebreaking process implemented in the algorithm calculates the minimal RMSD between two molecules.

## Usage notes

### Usage
 `Canonize.py [-h] [-o OUTPUT] file`

to obtain canonized indices and stereochemical tags for atoms in a molecule.

 `CanonizedRMSD.py [-h] [-s] [-m] [-i] [-q] [-r] file1 file2`

to calculate the RMSD of two molecules after canonizing them.
<br>

### Supported file types
 `.mol | .sdf | .rxn | .mol2 | .ml2 | .pdb`


### Optional arguments:

<pre>

Canonize.py file
    -o OUTPUT, --output OUTPUT    output file name. When not specified, the output will be printed to the screen.

CanonizedRMSD.py file1 file2

    -s, --save save intermediate results

    -m, --mapping output atom mapping relationship with two molecules

    -i, --ignore_isomerism ignore geometric and stereometric isomerism when canonizing

    -a, --no_alignment do not apply molecule alignment by Kabsch algorithm or QCP algorithm when calculating RMSD

    -q, --use QCP algorithm instead of Kabsch algorithm

    -r, --remove H atoms

    -at, --arbitrary_tiebreaking  apply an arbitrary tiebreaking, namely not performing branching tiebreaking
</pre>