from rdkit import Chem
from typing import Literal, Union  

from CanonizedRMSD.utils import prepare_molecule
from CanonizedRMSD.core import get_canonized_mapping, calc_canonical_rmsd
from CanonizedRMSD.data import CanonizedMapping, CanonizedRMSDResult


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

    def canonize(self, input_file: Union[str, Chem.rdchem.Mol]) -> CanonizedMapping:
        molecule, _ = prepare_molecule(input_file, 
                                       removeH=False, 
                                       aromatize=self.aromatize, 
                                       assign_rdkit_stereo=self.rdkit_stereo)
        molecule_noH = Chem.RemoveHs(molecule)
        mapping = get_canonized_mapping(molecule_noH, stereo=self.rdkit_stereo)
        return mapping

class RMSDCalc:
    '''
    A calculator for canonized and symmetry-corrected RMSD calculation between two molecules.
    '''
    def __init__(self, ignore_isomerism: bool=False,
                 identity_check: bool=False,
                 remove_Hs: bool=False,
                 rmsd_algorithm: Literal["Kabsch", "QCP"]="Kabsch",
                 symmetry_correction: bool=True) -> None:
        self.ignore_isomerism = ignore_isomerism
        self.identity_check = identity_check
        self.remove_Hs = remove_Hs
        self.rmsd_algorithm = rmsd_algorithm
        self.symmetry_correction = symmetry_correction

    def run(self, input_file_1: Union[str, Chem.rdchem.Mol],
                  input_file_2: Union[str, Chem.rdchem.Mol],
                  no_alignment: bool=False) -> CanonizedRMSDResult:
        mol_A, non_H_idx_A = prepare_molecule(input_file_1, removeH=self.remove_Hs)
        mol_B, non_H_idx_B = prepare_molecule(input_file_2, removeH=self.remove_Hs)
        canonized_result = calc_canonical_rmsd(mol_A, 
                                               mol_B, 
                                               no_isomerism=self.ignore_isomerism,
                                               identity_check=self.identity_check,
                                               no_alignment=no_alignment,
                                               algorithm=self.rmsd_algorithm,
                                               tiebreaking=self.symmetry_correction)
        # specify the heavy atom index remappings to make the index work properly when remove_Hs=True
        if self.remove_Hs:
            canonized_result.file1_mapping.non_h_idx_mapping = non_H_idx_A
            canonized_result.file2_mapping.non_h_idx_mapping = non_H_idx_B
            canonized_result.file1_mapping._remap_non_H_idx()
            canonized_result.file2_mapping._remap_non_H_idx()
        return canonized_result