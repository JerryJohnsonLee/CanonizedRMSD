#!/usr/bin/env python
import sys
import argparse
import time

import numpy as np

try:
    from rdkit import Chem
    import formatting
    import main
except:
    print("RDKit module not found in your Python! Please make sure RDKit is correctly installed.")
    print("Please go to http://www.rdkit.org/docs/Install.html for more details about installation.")
    sys.exit()

stereochemical_tags = ['/', 's', 'r', 'S', 'R', 'E', 'Z', ':', '?']

def canonize(file_path, aromatize=False, stereo=False):
    fileState = main.CheckValidity(file_path)
    molecule, _ = formatting.Read(file_path, 0, fileState, False,
                aromatize=aromatize, assignRDKitStereo=stereo)
    molecule_noH = Chem.RemoveHs(molecule)
    content, _ = main.CanonizedSequenceRetriever(molecule_noH, stereo=stereo)
    return sorted(content,key=lambda x:x['original'])

if __name__=="__main__":
    global fileState
    DEBUG=False
    if not DEBUG:
        parser=argparse.ArgumentParser( \
        description="to obtain canonized indices and stereochemical tags for atoms in a molecule. \
        \n\nsupported file types:\n   .mol | .sdf | .rxn | .mol2 | .ml2  | .pdb\n", \
        formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument("file")
        parser.add_argument("-o", "--output", default=None,\
             help="output file name. When not specified, the output will be printed to the screen.")
        parser.add_argument("-r", "--aromatize", action="store_true", 
             help="aromatize the molecule (which may change its structure) before canonization. Default to False.")
        parser.add_argument("-s", "--stereo", action="store_true", 
             help="assign stereo according to RDKit tags prior to canonization, which may help with dependent stereochemistry. Default to False.")

        args = parser.parse_args()
        canonization = canonize(args.file, args.aromatize, args.stereo)

        output_results = ["Original Index     Canonized Index     Atom Type   Stereochemical Tag"]
        output_results.append("-" * 80)
        for item in canonization:
            output_results.append("{:>14}     {:>15}     {:>9}   {:>18}".format(\
                item['original'] + 1, item['canonized'] + 1,\
                     item['item'].source.GetSymbol(), stereochemical_tags[item['item'].stereochemistry]))

        print("\n".join(output_results))
        if args.output is not None:
            with open(args.output, 'w') as f:
                f.write("\n".join(output_results))

