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

stereochemical_tags = ['/', 's', 'r', 'S', 'R', 'E', 'Z']

def canonize(file_path):
    fileState = main.CheckValidity(file_path)
    molecule, _ = formatting.Read(file_path, 0, fileState, False)
    content, _ = main.CanonizedSequenceRetriever(molecule)
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

        args = parser.parse_args()
        canonization = canonize(args.file)

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

