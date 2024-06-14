'''Define the various data structures used in the project, and their special methods'''
from rdkit import Chem
from typing import List
from dataclasses import dataclass
import pandas as pd

import numpy as np

R_TYPE = 4
S_TYPE = 3
Z_TYPE = 6
E_TYPE = 5
r_TYPE = 2
s_TYPE = 1
AROMATIC = 7
UNKNOWN = 8
DOUBLE_BOND = Chem.rdchem.BondType.DOUBLE
AROMATIC_BOND = Chem.rdchem.BondType.AROMATIC
stereochemical_tags = ['/', 's', 'r', 'S', 'R', 'E', 'Z', ':', '?']

@dataclass
class RMSDResult:
    rmsd: float # The RMSD value
    transition: np.ndarray = None  # The translation vector
    rotation: np.ndarray = None  # The rotation matrix
    transformed: np.ndarray = None  # The transformed coordinates

@dataclass
class CanonizedMapping:
    lst: List

    def __post_init__(self) -> None:
        self._sort()

    def _sort(self) -> None:
        '''Sort the canonized mapping in place.'''
        self.lst = sorted(self.lst, key=lambda x: x['original'])
        return self
    
    def __str__(self) -> str:
        '''get the string representation of the canonized mapping.'''
        output_results = ["Original Index     Canonized Index     Atom Type   Stereochemical Tag"]
        output_results.append("-" * 80)
        for item in self.lst:
            output_results.append("{:>14}     {:>15}     {:>9}   {:>18}".format(\
                item['original'] + 1, item['canonized'] + 1,\
                item['item'].source.GetSymbol(), stereochemical_tags[item['item'].stereochemistry]))
            
        result = "\n".join(output_results)
        return result
    
    def dataframe(self) -> pd.DataFrame:
        '''Return the canonized mapping as a pandas DataFrame.'''
        return pd.DataFrame(self.lst)
    
    def save(self, path: str) -> None:
        '''Save the canonized mapping to a file.'''
        with open(path, 'w') as f:
            f.write(str(self))
