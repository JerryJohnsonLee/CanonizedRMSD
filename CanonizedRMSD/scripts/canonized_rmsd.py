"""
Entry point script for calculating canonized RMSDs between two molecules

Author: Jie Li
Date created: Jun 25, 2024
"""

import time
import argparse
from CanonizedRMSD.main import RMSDCalc

def parse_args():
    '''Prepare argument parser'''
    parser = argparse.ArgumentParser( \
        description="to calculate the RMSD of two molecules after canonizing them. \
    \n\nsupported file types:\n   .mol | .sdf | .rxn | .mol2 | .ml2  | .pdb\n", \
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("file1")
    parser.add_argument("file2")
    parser.add_argument("-s", "--save", action="store_true", help="save intermediate results")
    parser.add_argument("-m", "--mapping", action="store_true", help="output atom mapping relationship with two molecules")
    parser.add_argument('-i', "--ignore_isomerism", action="store_true", help="ignore geometric and stereometric isomerism when canonizing")
    parser.add_argument('-ic', "--identity_check", action="store_true", help="perform molecule identity check based on non-hydrogen graph before canonizing")
    parser.add_argument('-na', "--no_alignment", action="store_true", help="do not apply molecule alignment Kabsch or QCP algorithm when calculating RMSD")
    parser.add_argument('-alg', "--algorithm", choices=["Kabsch", "QCP"], default="Kabsch", help="algorithm to calculate RMSD, default to Kabsch")
    parser.add_argument('-r', "--removeHs", action="store_true", help="remove H atoms")
    parser.add_argument('-nsc', "--no_symmetry_correction", action="store_true", help="not performing symmetry correction (branching tiebreaking) during RMSD calculation")
    parser.add_argument('-t', '--time', action='store_true', help='print total calculation time')

    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    calculator = RMSDCalc(args.ignore_isomerism, args.identity_check, args.removeHs, args.algorithm, not args.no_symmetry_correction)
    if args.time:
        t0 = time.time()
    results = calculator.run(args.file1, args.file2, args.no_alignment)

    if args.time:
        print("Calculation finished in %f s." % (time.time() - t0))
    print("Canonized RMSD:", results.rmsd)

    if args.mapping:
        print("=" * 40)
        print(results.mapping_repr)

if __name__ == "__main__":
    main()