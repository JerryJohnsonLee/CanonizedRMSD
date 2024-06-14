'''Utility functions to do format conversions, structure conversions, etc.'''
from io import StringIO
from typing import List
import numpy as np


def remove_Hs_from_mol_block(input_mol_block: str) -> str:
    '''
    Remove hydrogen atoms from a molecule represented as block string.

    Args:
        input_mol_block: str, the content of the molecule block

    Returns:
        string, the content of the file without hydrogen atoms
        sequence: list, indices of the heavy atoms in the original molecule
    '''
    s_file = StringIO()
    s_file.write(input_mol_block)
    s_file.seek(0)

    content = []
    content.append("\n     RDKit          \n\n")

    s_file.readline()
    s_file.readline()
    s_file.readline()
    l = s_file.readline()
    
    n_atoms = int(l[:3])
    atoms = []
    n_bonds = int(l[3:6])
    sequence = []
    index = 1
    for _ in range(0, n_atoms):
        line = s_file.readline()
        atoms.append(line)
        if 'H  ' != line[31:34]:
            content.append(line)
            sequence.append(index)
            index += 1
        else:
            sequence.append(0)
    bond_counts = 0
    for _ in range(0, n_bonds):
        l = s_file.readline()
        n1 = int(l[:3])
        n2 = int(l[3:6])
        if sequence[n1 - 1] * sequence[n2 - 1] != 0:
            content.append("%3s%3s" % (sequence[n1 - 1],
                                       sequence[n2 - 1]) + l[6:12] + '\n')
            bond_counts += 1
    atom_counts = n_atoms - sequence.count(0)
    content.insert(1, "%3d%3d" % (atom_counts, bond_counts) + \
                   "  0  0  0  0  0  0  0  0999 V2000\n")
    content.append("M  END\n$$$$\n")
    string = ''.join(content)
    return (string, sequence)

def check_elements(e1: List[str], e2: List[str]) -> bool:
    '''
    Check if two lists of elements are identical.

    Args:
        e1: list, a list of elements provided as strings
        e2: list, a list of elements provided as strings
    '''
    if len(e1) != len(e2):
        return False
    else:
        for i in range(0, len(e1)):
            if e1[i] != e2[i]:
                return False
        return True

def get_coords_elements(mol_block: str):
    '''
    Get the coordinates as numpy array and elements from a molecule block.
    '''
    lines = []
    points = []
    elements = []
    for line in mol_block.split('\n'):
        lines.append(line)

    n_atoms = int(lines[3][:3])
    for n in range(4, 4 + n_atoms):
        l = lines[n]
        points.append([float(l[:10]), float(l[10:20]), float(l[20:30])])
        elements.append(l[31:34].strip())
    return np.matrix(points), elements