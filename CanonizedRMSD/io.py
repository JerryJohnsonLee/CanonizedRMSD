'''
Helper function module for inputs and outputs of certain molecular formats
'''

from typing import Optional, List
from rdkit import Chem

def get_non_hydrogen_indices(molecule: Chem.Mol) -> List[int]:
    # get the non-H atom indices and form a list
    sequence = []
    for atom in molecule.GetAtoms():
        if atom.GetSymbol() != "H":
            sequence.append(atom.GetIdx())
    return sequence


def convert_from_gaussian_to_rdkit(input_file: str, output_file: Optional[str]=None) -> str:
    '''
    Do format conversion from Gaussian to a format readable by rdkit. 
    The RDKit format should have "  0999 V2000" in the 4th line.

    Args:
        input_file: str, the path to the input file
        output_file: str, the path to the output file to be saved if it is not None

    Returns:
        string, the content of the converted file
    ''' 
    try:  
        s_file = open(input_file, 'r')
    except:
        print("Cannot open file %s, please check!" % input_file)
        import sys
        sys.exit()
    content = []
    content.append("\n     RDKit          \n\n")
    s_file.readline()
    s_file.readline()
    s_file.readline()
    l = s_file.readline()

    if(len(l) == 40):    # already rdkit type
        s_file.seek(0)
        return s_file.read()
    content.append(l[:27] + "  0999 V2000\n")

    content.append(s_file.read())
    content.append("M  END\n$$$$\n")
    
    string = ''.join(content)

    # when provided output file path, write the content to the file
    if  output_file is not None:
        t_file = open(output_file,'w')
        t_file.write(string)
        t_file.close()
    s_file.close()
    return string

def _read_from_mol(file: str, output_filepath: Optional[str]=None) -> Chem.Mol:
    converted_mol_block = convert_from_gaussian_to_rdkit(file, output_filepath)
    molA = Chem.MolFromMolBlock(converted_mol_block, sanitize=False)
    return molA

def _read_from_mol2(file: str) -> Chem.Mol:
    molA = Chem.MolFromMol2File(file, sanitize=False)
    return molA

def _read_from_pdb(file: str) -> Chem.Mol:
    molA = Chem.MolFromPDBFile(file, sanitize=False)
    return molA

mol_readers = {
    "mol": _read_from_mol,
    "sdf": _read_from_mol,
    "mol2": _read_from_mol2,
    "pdb": _read_from_pdb
}

def _parse_extension_from_filename(file: str) -> str:
    return file.split(".")[-1]

def read(file: str, removeH: bool, aromatize: bool=False, assign_rdkit_stereo: bool=False):
    '''
    Main function for reading a molecule from a file.

    Args:
        file: str, the path to the input file
        removeH: bool, whether to remove H atoms
        aromatize: bool, whether to aromatize the molecule
        assign_rdkit_stereo: bool, whether to assign stereochemistry tags using RDKit 
            Chem.AssignStereochemistryFrom3D method

    Returns:
        mol: RDKit Mol object, the molecule read from the file
        non_H_indices: list, the indices of the non-H atoms in the molecule
    '''
    extension = _parse_extension_from_filename(file)
    if extension in mol_readers:
        mol = mol_readers[extension](file)
    else:
        try:
            mol = _read_from_mol(file)
        except:
            try:
                mol = _read_from_mol2(file)
            except:
                print("Unsupported file: %s" % file)
                import sys
                sys.exit()

    # remove Hs when required
    if removeH:
        mol=Chem.RemoveHs(mol)

    # partial sanitization to make sure necessary molecular properties are calculated
    mol.UpdatePropertyCache(strict=False)
    sanitizeOpts=Chem.SanitizeFlags.SANITIZE_FINDRADICALS| \
                 Chem.SanitizeFlags.SANITIZE_KEKULIZE| \
                 Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION| \
                 Chem.SanitizeFlags.SANITIZE_SYMMRINGS| \
                 Chem.SanitizeFlags.SANITIZE_ADJUSTHS
    if aromatize:
        sanitizeOpts = sanitizeOpts|Chem.SanitizeFlags.SANITIZE_SETAROMATICITY| \
                                    Chem.SanitizeFlags.SANITIZE_SETCONJUGATION

    sanitization_result = Chem.SanitizeMol(mol,sanitizeOps=sanitizeOpts,
                                                catchErrors=True)
    if sanitization_result != Chem.SanitizeFlags.SANITIZE_NONE:
        print("Cannot sanitize file: %s" % file)
        import sys
        sys.exit()

    # assign rdkit stereochemistry tags to be used for initializing idCode (when necessary)
    if assign_rdkit_stereo:
        Chem.AssignStereochemistryFrom3D(mol)

    non_H_indices = get_non_hydrogen_indices(mol)
    return mol, non_H_indices