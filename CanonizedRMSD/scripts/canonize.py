"""
Entry point script for canonizing one molecule

Author: Jie Li
Date created: Jun 14, 2024
"""

import argparse
from CanonizedRMSD.main import Canonizer

def parse_args():
    '''Prepare argument parser'''
    parser = argparse.ArgumentParser( \
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
    return args

def main():
    args = parse_args()
    canonizer = Canonizer(args.aromatize, args.stereo)
    results = canonizer.canonize(args.file)

    print(results)
    if args.output is not None:
        results.save(args.output)

if __name__ == "__main__":
    main()