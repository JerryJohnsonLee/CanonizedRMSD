import pytest
import os

from CanonizedRMSD import RMSDCalc
from rdkit import Chem
import random

calc = RMSDCalc()
molecules = os.listdir(os.path.join(pytest.TEST_FOLDER, "VS_compounds"))

def random_renumber(input_file_path):
    mol = Chem.MolFromMolFile(input_file_path, removeHs=False)
    atom_num = mol.GetNumAtoms()
    new_order = list(range(atom_num))
    random.shuffle(new_order)
    new_mol = Chem.RenumberAtoms(mol, newOrder=new_order)
    return mol, new_mol

@pytest.mark.parametrize("name", molecules)
def test_vs_compounds(name):
    mol, new_mol = random_renumber(f"{pytest.TEST_FOLDER}/VS_compounds/{name}")
    result = calc.run(mol, new_mol)
    assert pytest.approx(result.rmsd, abs=1e-6) == 0
