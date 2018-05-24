# augmented_Kabsch
The augmented Kabsch algorithm with canonization


usage: CanonizedRMSD.py [-h] [-s] [-m] [-i] file1 file2

to calculate the RMSD of two molecules after canonizing them.

supported file types:
   .mol | .sdf | .rxn | .mol2 | .ml2

positional arguments:
  file1
  file2

optional arguments:

  -h, --help            show this help message and exit
  
  -s, --save            save intermediate results
  
  -m, --mapping         output atom mapping relationship with two molecules
  
  -i, --ignore_isomerism
                        ignore geometric and stereometric isomerism when canonizing
                        
  -a, --no_alignment    do not apply molecule alignment by Kabsch algorithm when calculating RMSD

