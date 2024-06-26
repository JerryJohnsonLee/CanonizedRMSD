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

    def __gt__(self, other: 'RMSDResult') -> bool:
        '''Compare two RMSD results based on the RMSD value.'''
        return self.rmsd > other.rmsd
    
    def __lt__(self, other: 'RMSDResult') -> bool:
        '''Compare two RMSD results based on the RMSD value.'''
        return self.rmsd < other.rmsd

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

    @classmethod
    def from_atom_list(cls, atom_list: List[atom]) -> 'CanonizedMapping':
        '''Create a CanonizedMapping object from a list of Atoms.'''
        return cls([{"original": atom.originalIndex,
                     "canonized": atom.currentIndex, 
                     "item": atom} for atom in atom_list])

    @property
    def original_to_canonized_mapping(self) -> List:
        '''Get the mapping from original index to canonized index.'''
        return [i["canonized"] for i in sorted(self.lst, key=lambda p:p["original"])]
    
    @property
    def canonized_to_original_mapping(self) -> List:
        '''Get the mapping from canonized index to original index.'''
        return [i["original"] for i in sorted(self.lst, key=lambda p:p["canonized"])]

@dataclass
class CanonizedRMSDResult(RMSDResult):
    file1_mapping: CanonizedMapping = None
    file2_mapping: CanonizedMapping = None

    @classmethod
    def from_rmsd_result(cls, rmsd_result: RMSDResult):
        '''Update the RMSD result with the values from the RMSD result.'''
        return cls(rmsd=rmsd_result.rmsd, 
                   transition=rmsd_result.transition, 
                   rotation=rmsd_result.rotation, 
                   transformed=rmsd_result.transformed)

class atom:
    def __init__(self,a=None,OriginalIndex=0,CopyFrom=None,addRDKitStereo=False,minAtomRingSize=0,coord=None):
        if CopyFrom==None:
            self.source=a
            self.degree=a.GetDegree()
            self.atomicNumber=a.GetAtomicNum()
            self.attatchedHs=a.GetTotalNumHs(includeNeighbors=True)
            self.charge=a.GetNumRadicalElectrons()
            self.isotope=a.GetIsotope()
            self.minAtomRingSize=minAtomRingSize
            self.stereochemistry=0
            if addRDKitStereo:
                self.rdkitStereoTag=int(a.GetChiralTag())
            else:
                self.rdkitStereoTag=0
            self.bonds=a.GetBonds()
            self.bondTypes=[bond.GetBondType() for bond in a.GetBonds()]
            self.doubleBondConnectAtom=None
            self.neighbors=()
            self.neighborIndexs=[]
            self.neighborElements=()
            self.originalIndex=OriginalIndex
            self.currentIndex=0
            self.isComplete=False
            self.coordinate=coord
        else:
            self.source=CopyFrom.source
            self.neighbors=()
            self.neighborIndexs=[]
            self.originalIndex=CopyFrom.originalIndex
            self.currentIndex=CopyFrom.currentIndex
            self.atomicNumber=CopyFrom.atomicNumber
            self.bonds=CopyFrom.bonds
            self.isComplete=CopyFrom.isComplete
            self.stereochemistry=CopyFrom.stereochemistry
            self.doubleBondConnectAtom=CopyFrom.doubleBondConnectAtom
            self.idcode=CopyFrom.idcode
            self.minAtomRingSize=CopyFrom.minAtomRingSize
            self.coordinate=CopyFrom.coordinate

    def assignIdCode(self, loose_criterion=False, no_H=False):
        # new criterion:
        #    ***                ***                  *       *       *        *           *
        # atomicNum    largestNeighborElement   connection  n_Hs  isotope   charge  rdkitStereoTag
        if loose_criterion:
            connection=self.degree
        else:
            connection=self.degree+self.bondTypes.count(Chem.rdchem.BondType.DOUBLE) \
            +self.bondTypes.count(Chem.rdchem.BondType.TRIPLE)*2 \
            +self.bondTypes.count(Chem.rdchem.BondType.AROMATIC)*0.5
        if len(self.neighbors)==0:
            largestNeighborElement=0  # i.e. single ions or single atoms
        else:
            largestNeighborElement=max([item.atomicNumber for item in self.neighbors])
        idCode=0
        idCode+=self.atomicNumber*1e8
        idCode+=largestNeighborElement*1e5
        idCode+=connection*1e4
        if not no_H:
            idCode+=self.attatchedHs*1e3
        idCode+=self.isotope*100
        if not loose_criterion:
            idCode+=(self.charge+4)*10   # to avoid problem with negative charges
        idCode+=self.rdkitStereoTag
        self.idcode=int(idCode)

    def completeNeighborElements(self):
        self.neighborElements=tuple(sorted([neighborAtom.source.GetAtomicNum() \
     for neighborAtom in self.neighbors],reverse=True))   # a tuple is used to make it hashable

    def updateNeighborIndexs(self):
        self.neighborIndexs=[item.currentIndex for item in self.neighbors]
        self.neighborIndexs.sort(reverse=True)

    def decideDoubleBondConnectedAtom(self):
        if self.bondTypes.count(DOUBLE_BOND) == 1:
        # if self.atomicNumber==6 and DOUBLE_BOND in self.bondTypes:
            doubleBond=[bond for bond in self.bonds if bond.GetBondType()==DOUBLE_BOND][0]
            doubleBondAtom=[item for item in self.neighbors \
                if item.source.GetIdx() in (doubleBond.GetBeginAtomIdx(),doubleBond.GetEndAtomIdx())][0]
            self.doubleBondConnectAtom=doubleBondAtom

