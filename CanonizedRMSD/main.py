from rdkit import Chem

from CanonizedRMSD.io import read
from CanonizedRMSD.core import get_canonized_mapping
from CanonizedRMSD.data import CanonizedMapping


class Canonizer:
    '''
    The main class for canonizing molecules. It is initialized with two parameters:
    aromatize: whether to aromatize the input molecules, defaults to False
    rdkit_stereo: whether to assign stereochemistry tags using RDKit Chem.AssignStereochemistryFrom3D method,
        This might provide additional information when initializing the ID code
        for the canonization algorithm. Defaults to False
    '''
    def __init__(self, aromatize: bool=False, rdkit_stereo: bool=False) -> None:
        self.aromatize = aromatize
        self.rdkit_stereo = rdkit_stereo

    def canonize(self, input_file: str) -> CanonizedMapping:
        molecule, _ = read(input_file, removeH=False,
                            aromatize=self.aromatize, assign_rdkit_stereo=self.rdkit_stereo)
        molecule_noH = Chem.RemoveHs(molecule)
        mapping = get_canonized_mapping(molecule_noH, stereo=self.rdkit_stereo)
        return mapping

class RMSDCalc:
    def __init__(self) -> None:
        pass