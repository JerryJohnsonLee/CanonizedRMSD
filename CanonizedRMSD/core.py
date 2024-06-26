from rdkit import Chem
import numpy as np
import formatting
import os
import sys
from typing import Literal, Tuple

from CanonizedRMSD.data import *
from CanonizedRMSD.utils import check_identity, check_elements
from CanonizedRMSD.rmsd import RMSD_CALCULATORS, calc_rmsd_with_mapping


def get_canonized_mapping(
        mol: Chem.rdchem.Mol,
        no_isomerism: bool=False,
        unbroken_molecule_container=[],
        no_Hs=False,
        stereo=False) -> CanonizedMapping:
    if len(unbroken_molecule_container) == 0:
        assert mol is not None, "Molecule is None"
        unbroken_molecule = canonize_molecule(mol,no_isomerism,no_Hs,stereo=stereo)
        unbroken_molecule_container.append(unbroken_molecule)
    # if serial:
    #     return branching_tiebreaking(mol,unbroken_molecule,ma,ea,no_alignment,qcp)

    collection = deepcopy(unbroken_molecule)
    ordinary_tiebreaking(collection)
    mapping = CanonizedMapping.from_atom_list(collection)
    return mapping
    
def check_noH_identity(mol_A: Chem.rdchem.Mol,
                       mol_B: Chem.rdchem.Mol,
                       no_isomerism: bool=False) -> bool:
    mol_A_noH = Chem.RemoveHs(mol_A)
    mol_B_noH = Chem.RemoveHs(mol_B)
    mapping_A_noH = get_canonized_mapping(mol_A_noH, no_isomerism=no_isomerism, no_Hs=True)
    mapping_B_noH = get_canonized_mapping(mol_B_noH, no_isomerism=no_isomerism, no_Hs=True)
    identity = check_identity(mapping_A_noH, mapping_B_noH)
    return identity

def calc_canonical_rmsd(mol_A: Chem.rdchem.Mol,
                        mol_B: Chem.rdchem.Mol, 
                        no_isomerism: bool=False, 
                        identity_check: bool=False, 
                        no_alignment: bool=False, 
                        algorithm: Literal["Kabsch", "QCP"]="Kabsch", 
                        tiebreaking: bool=True) -> CanonizedRMSDResult:
    mapping_A = get_canonized_mapping(mol_A, no_isomerism=no_isomerism, unbroken_molecule_container=[])
    # prepare a container for the unbroken molecule so that it can be used later
    unbroken_molecule_container = []
    temp_mapping_B = get_canonized_mapping(mol_B, no_isomerism=no_isomerism,
                                      unbroken_molecule_container=unbroken_molecule_container)

    if identity_check:
        if not check_noH_identity(mol_A, mol_B, no_isomerism):
            raise RuntimeError("The two input molecules are not identical based on non-hydrogen molecule graph!")
    canonical_mol_A = Chem.RenumberAtoms(mol_A, newOrder=mapping_A.canonized_to_original_mapping)
    coord_mat_A = canonical_mol_A.GetConformer().GetPositions()
    elements_A = [atm.GetSymbol() for atm in canonical_mol_A.GetAtoms()]

    if tiebreaking:
        min_RMSD_mapping_B, rmsd_result = branching_tiebreaking(mol_B, unbroken_molecule_container[0], coord_mat_A, elements_A, 
                                                     no_alignment=no_alignment, algorithm=algorithm)
    else:
        min_RMSD_mapping_B = temp_mapping_B
        rmsd_result = calc_rmsd_with_mapping(coord_mat_A, elements_A, mol_B, min_RMSD_mapping_B, no_alignment=no_alignment, algorithm=algorithm)
    canonized_result = CanonizedRMSDResult.from_rmsd_result(rmsd_result)
    canonized_result.file1_mapping = mapping_A
    canonized_result.file2_mapping = min_RMSD_mapping_B
    return canonized_result


def branching_tiebreaking(
    mol_B: Chem.rdchem.Mol,
    unbroken_molecule: List[atom],
    coord_mat_A: np.ndarray,
    elements_A: List,
    no_alignment: bool=False,
    algorithm: Literal["Kabsch", "QCP"]="Kabsch") -> Tuple[CanonizedMapping, RMSDResult]:
    unbroken_idxs = []
    unbroken_idx_2_original = {}
    for item in unbroken_molecule:
        if not item.isComplete:
            unbroken_idxs.append(item.currentIndex)
            AddItem(unbroken_idx_2_original, item.currentIndex, item.originalIndex)
    # when there is no unbroken index, simply calculate RMSD and return
    if len(unbroken_idxs) == 0:
        mapping = CanonizedMapping.from_atom_list(unbroken_molecule)
        rmsd_result = calc_rmsd_with_mapping(coord_mat_A, elements_A, mol_B, mapping, no_alignment=no_alignment, algorithm=algorithm)
        return mapping, rmsd_result

    selected_idx = select_best_idx_for_tiebreaking(unbroken_idxs)
    selected_number = len(unbroken_idx_2_original[selected_idx])
    items_with_chosen_idx = [item for item in unbroken_molecule if item.currentIndex == selected_idx]

    # enumerate all possibilities of selecting the first atom to assign complete tag
    min_rmsd = RMSDResult(99999)
    min_rmsd_mapping = None
    min_rmsd_atom_lst = None

    for i in range(selected_number):
        molecule_tmp = deepcopy(unbroken_molecule)
        items_with_chosen_idx = [item for item in molecule_tmp if item.currentIndex == selected_idx]
        for idx in range(selected_number):
            if idx == i:
                items_with_chosen_idx[idx].isComplete = True
            else:
                # for the rest, increase their index by 1
                items_with_chosen_idx[idx].currentIndex += 1
        
        working_set = set()
        for item in molecule_tmp:
            if not item.isComplete:
                item.updateNeighborIndexs()
            working_set.add(item)
        
        Refine(working_set)
        tmp_mapping = CanonizedMapping.from_atom_list(molecule_tmp)
        rmsd_result = calc_rmsd_with_mapping(coord_mat_A, elements_A, mol_B, tmp_mapping, 
                                             no_alignment=no_alignment, algorithm=algorithm, completed_atoms_only=True)

        if rmsd_result < min_rmsd:
            min_rmsd = rmsd_result
            min_rmsd_mapping = tmp_mapping
            min_rmsd_atom_lst = molecule_tmp

    
    # finally, decide whether the tiebreaking is complete (so that the results can be returned back to the parent function),
    # or recursive treatment is needed
    complete = True
    for item in molecule_tmp:
        if not item.isComplete:
            complete=False
            break

    if complete:
        return min_rmsd_mapping, min_rmsd
    else:
        return branching_tiebreaking(mol_B, min_rmsd_atom_lst, coord_mat_A, elements_A,
                                     no_alignment=no_alignment, algorithm=algorithm)

def deepcopy(worklist):
    newlist = []
    for item in worklist:
        newlist.append(atom(CopyFrom=item))
    CompleteNeighbor(newlist)
    return newlist


def CompleteNeighbor(worklist):
    # Form neighbor network in worklist
    for item in worklist:
        item.neighbors=tuple(p for p in worklist if (p.source.GetIdx() in [q.GetIdx() for q in item.source.GetNeighbors()]))
        item.updateNeighborIndexs()
        neighborMapDict={neighbor.originalIndex:neighbor for neighbor in item.neighbors}
        item.bondConnectedAtoms=[neighborMapDict[bond.GetOtherAtomIdx(item.originalIndex)] for bond in item.bonds]
    

def AddItem(dictionary,key,value):
    if dictionary.__contains__(key):
        dictionary[key].append(value)
    else:
        dictionary[key]=[value]

def SelectBestIndex(workset):
    # Select the partition with most instances as current partition
    currentIndexs=[item.currentIndex for item in workset]
    indexNumbers=[]
    for number in set(currentIndexs):
        indexNumbers.append({"index":number,"count":currentIndexs.count(number)})
    selectedItem=sorted(indexNumbers,key=lambda x:(-x["count"],-x["index"]))[0]
    return selectedItem["index"]



def Refine(workset):
    dictionaryForIndexs={}
    for item in workset:
        AddItem(dictionaryForIndexs,item.currentIndex,item)
    currentIndexMap=constructCurrentIndexMap(workset)
    while len(workset)>0:        
        currentIndex=SelectBestIndex(workset)
        currentPartition=[items for items in dictionaryForIndexs[currentIndex] if items in workset]
        currentPartition.sort(key=lambda i:i.neighborIndexs)
        changedIndex=False
        lastIndexCollection=currentPartition[0].neighborIndexs
        i=0
        for (index,item) in enumerate(currentPartition):
            if item.neighborIndexs!=lastIndexCollection:
                lastIndexCollection=item.neighborIndexs
                i=index
            if i!=0:
                dictionaryForIndexs[currentIndex].remove(item)
                item.currentIndex=currentIndex+i
                AddItem(dictionaryForIndexs,currentIndex+i,item)
                changedIndex=True
        # Decide whether there are completed atoms
        for item in currentPartition:
            if len(dictionaryForIndexs[item.currentIndex])==1 and not item.isComplete:
                item.isComplete=True
                #workset.remove(item)
                # print(item.originalIndex,item.currentIndex,"becomes True in A") #
        # When no element in the current partition needs index change
        #if not changedIndex:
            # Remove these elements from workset
        for item in currentPartition:
            workset.remove(item)
        # When changes in current partition indexs
        #else:
        if changedIndex:
            # Update neighbor indexs of changed atoms' neighbors and recheck them
            for item in currentPartition:
                for neighbor in item.neighbors:
                    neighbor.updateNeighborIndexs()
                    if dictionaryForIndexs.__contains__(neighbor.currentIndex):
                        for item in dictionaryForIndexs[neighbor.currentIndex]:
                            if not item.isComplete:
                                addItemWithCurrentIndexMap(item,workset,currentIndexMap)

def select_best_idx_for_tiebreaking(unbrokenIndexs):
    indexNumbers=[]
    for index in set(unbrokenIndexs):
        indexNumbers.append({'index':index,
                                'number':unbrokenIndexs.count(index)})
    # Select the index with smallest minimum topological distance, most occuring number while being as large as possible
    selectedItem=sorted(indexNumbers,key=lambda p:(-p['number'],-p['index']))[0]
    selectedIndex=selectedItem['index']
    return selectedIndex



def ordinary_tiebreaking(worklist):
    unbrokenIndexs = []
    for item in worklist:
        if not item.isComplete:
            unbrokenIndexs.append(item.currentIndex)
    if len(unbrokenIndexs) != 0:
        selectedIndex = select_best_idx_for_tiebreaking(unbrokenIndexs)

        # choose one item to assign complete tag. For all the rest ones with the same index,
        # assign +1 index. Then do refinement
        itemsWithChosenIndex=[item for item in worklist if item.currentIndex==selectedIndex]
        itemsWithChosenIndex[0].isComplete=True
        for idx in range(1, len(itemsWithChosenIndex)):
            itemsWithChosenIndex[idx].currentIndex+=1

        newSet=set()
        for item in worklist:
            if not item.isComplete:
                item.updateNeighborIndexs()
            newSet.add(item)
        
        Refine(newSet)
        ordinary_tiebreaking(worklist)

def DecideStereo(centralAtom):
    neighborAtoms=centralAtom.neighbors
    if centralAtom.rdkitStereoTag==0:
        stereoTag2D=UNKNOWN  # we do not know 2D chiral tag from molecular graph
    else:
        # this is for sure a chiral center. Get the chiral tag from rdkit
        unsortedNeighborIndices=list(neighborAtom.currentIndex for neighborAtom in centralAtom.bondConnectedAtoms)
        sortedNeighborIndices=sorted(unsortedNeighborIndices)
        flip=False
        if len(unsortedNeighborIndices)==4:
            # all 4 connected atoms are explicit, then the last atom is the one behind the central atom
            behindAtomIdx=unsortedNeighborIndices[-1]
        else:
            behindAtomIdx=[idx for idx in centralAtom.neighborIndexs if idx not in unsortedNeighborIndices]
            if len(behindAtomIdx)==0:
                behindAtomIdx=0
            else:
                assert len(behindAtomIdx)==1
                behindAtomIdx=behindAtomIdx[0]
        if behindAtomIdx>min(sortedNeighborIndices):
            flip=not flip # flip if the atom behind the central atom is not the smallest index
        a,b,c=sortedNeighborIndices[:3]
        flip_configs=[(a,c,b),(b,a,c),(c,b,a)]
        if unsortedNeighborIndices in flip_configs:
            flip=not flip # decide whether the rdkit atom order needs to be flipped to comply with CIP
        if centralAtom.rdkitStereoTag==2:
            flip=not flip
        if flip:
            stereoTag2D=S_TYPE
        else:
            stereoTag2D=R_TYPE

    # get chiral tag from 3D coordinates
    if len(neighborAtoms)==4: # connected 4 different groups
        a,b,c,d=[np.array(neighborAtom.coordinate) for neighborAtom in sorted(neighborAtoms,key=lambda x:x.currentIndex,reverse=True)]
        # assign a,b,c,d to atomic coordinates according to their current index
        o=np.array(centralAtom.coordinate)
        ab=b-a
        ac=c-a
        od=d-o
        det=np.linalg.det(np.vstack([ab,ac,od]))
        if np.abs(det)<1e-5:
            stereoTag3D=UNKNOWN  # cannot infer chiral tag from 3D structure (i.e. it is actually only 2D)
        elif det>0: # ab x ac · od > 0 
            stereoTag3D=R_TYPE
        else:
            stereoTag3D=S_TYPE
    elif len(neighborAtoms)==3:
        a,b,c=[np.array(neighborAtom.coordinate) for neighborAtom in sorted(neighborAtoms,key=lambda x:x.currentIndex,reverse=True)]
        o=np.array(centralAtom.coordinate)
        ab=b-a
        ac=c-a
        oa=a-o
        det=np.linalg.det(np.vstack([ab,ac,oa]))
        if np.abs(det)<1e-5:
            stereoTag3D=UNKNOWN  # cannot infer chiral tag from 3D structure (i.e. it is actually only 2D)
        elif det>0:
            stereoTag3D=S_TYPE
        else:
            stereoTag3D=R_TYPE

    # finally, check consensus of 3d and 2d chiral tag, and output
    if stereoTag3D==UNKNOWN:
        return stereoTag2D
    elif stereoTag2D!=UNKNOWN and stereoTag3D!=stereoTag2D:
        print("Warning! Stereochemical tag for atom {} calculated from 3D coordinates is different from that encoded in the molecular graph. Stereochemical tag calculated from 3D coordinates is used.".format(centralAtom.originalIndex))
    return stereoTag3D

def ConvertedStereoTag(stereoTag,firstTime):
    if firstTime:
        return stereoTag
    else:
        if stereoTag==R_TYPE:
            return r_TYPE
        elif stereoTag==S_TYPE:
            return s_TYPE

def DecideDoubleBondGeometric(atom1,atom2):
    atom1Neighbors=set(atom1.neighbors)
    atom1Neighbors.remove(atom2)
    atom1BigNeighbor=sorted(atom1Neighbors,key=lambda p:p.currentIndex,reverse=True)[0]
    atom2Neighbors=set(atom2.neighbors)
    atom2Neighbors.remove(atom1)
    atom2BigNeighbor=sorted(atom2Neighbors,key=lambda p:p.currentIndex,reverse=True)[0]
    bondDirection1=np.array(atom1BigNeighbor.coordinate)-np.array(atom1.coordinate)
    bondDirection2=np.array(atom2BigNeighbor.coordinate)-np.array(atom2.coordinate)
    if bondDirection1.dot(bondDirection2)>0:
        return Z_TYPE
    else:
        return E_TYPE


def constructCurrentIndexMap(workset):
    indexMap={}
    for item in workset:
        if item.currentIndex in indexMap:
            indexMap[item.currentIndex].append(item)
        else:
            indexMap[item.currentIndex]=[item]
    return indexMap

def addItemWithCurrentIndexMap(item,workset,currentIndexMap):
    # make sure once an item is added back to workset, all items with the same current index are also added
    if item.currentIndex in currentIndexMap:
        itemsWithSameIndexAsNeighbor=currentIndexMap[item.currentIndex]
    else:
        itemsWithSameIndexAsNeighbor=[item]
    for content in itemsWithSameIndexAsNeighbor:
        workset.add(content)

def StereoAssigner(workset,firstTime=True):
    currentIndexMap=constructCurrentIndexMap(workset)
    while len(workset):
        currentIndex=SelectBestIndex(workset)
        for item in workset:
            if item.currentIndex==currentIndex:
                currentAtom=item
                break
        neighborCount=len(currentAtom.neighborIndexs)
        sameIndexAtoms=[item for item in workset if item.currentIndex==currentIndex] # find all atoms with same index
        updated=False
        if  currentAtom.stereochemistry==0: # check if the atom (carbon) has not been assigned stereo tag, start assigning
            if currentAtom.bondTypes.count(DOUBLE_BOND)==1: # one double bond (not in ring), not assigned geometric tag
                if currentAtom.minAtomRingSize not in [3, 4, 5, 6, 7]:  # only when the double bond is in a ring in size 3~7 will the algorithm not reporting its Z/E
                    if neighborCount>1 and len(set(currentAtom.neighborIndexs))==neighborCount:
                        doubleBondAtom=currentAtom.doubleBondConnectAtom
                        if (doubleBondAtom):
                            if len(doubleBondAtom.neighborIndexs)>1 and len(set(doubleBondAtom.neighborIndexs))==len(doubleBondAtom.neighborIndexs):
                                # both sides of the double bond have different groups
                                zCount=0
                                doubleBondAtoms=[]
                                for item in sameIndexAtoms:
                                    doubleBondAtom=item.doubleBondConnectAtom
                                    doubleBondAtoms.append(doubleBondAtom)
                                    geometricTag=DecideDoubleBondGeometric(item,doubleBondAtom)
                                    if geometricTag==Z_TYPE:
                                        zCount+=1
                                    item.stereochemistry=geometricTag
                                    doubleBondAtom.stereochemistry=geometricTag
                                updated=True
                                for atom1,atom2 in zip(sameIndexAtoms,doubleBondAtoms):
                                    if atom1 in workset:
                                        workset.remove(atom1)
                                    if atom2 in workset:
                                        workset.remove(atom2)
                                    
                                    if atom1.stereochemistry==Z_TYPE:
                                        # Z_TYPE atoms were updated their indices, while E_TYPE atoms keep their indices
                                        atom1.currentIndex+=len(sameIndexAtoms)-zCount
                                        atom2.currentIndex+=len(sameIndexAtoms)-zCount
                                        updateNeighborList=list(atom1.neighbors)
                                        updateNeighborList.extend(atom2.neighbors)
                                        for neighbor in updateNeighborList:
                                            neighbor.updateNeighborIndexs()
                                            if (not neighbor.isComplete) or neighbor.stereochemistry==0:
                                                addItemWithCurrentIndexMap(neighbor,workset,currentIndexMap)
                                if zCount:
                                    Refine(workset.copy())
            elif currentAtom.bondTypes.count(DOUBLE_BOND)==0 and currentAtom.bondTypes.count(AROMATIC_BOND)==0:  # no double bond or aromatic bond                  
                rCount=0        
                if neighborCount>2 and len(set(currentAtom.neighborIndexs))==neighborCount:
                # there is no same thing in current atom's neighbor list (implicit hydrogen included)                 
                    if currentAtom.atomicNumber==6 or (currentAtom.atomicNumber==7 and neighborCount==4): # carbon or tetrahedral nitrogen
                        for item in sameIndexAtoms:
                            stereoTag=DecideStereo(item)
                            if stereoTag==R_TYPE:
                                rCount+=1
                            item.stereochemistry=ConvertedStereoTag(stereoTag,firstTime)
                        updated=True
                        for item in sameIndexAtoms:
                            workset.remove(item)
                            if item.stereochemistry in (R_TYPE,r_TYPE):
                                item.currentIndex+=len(sameIndexAtoms)-rCount  # which is sCount
                            if not firstTime:
                                for neighbor in item.neighbors:
                                    neighbor.updateNeighborIndexs()
                                    if not neighbor.isComplete or neighbor.stereochemistry==0:
                                        addItemWithCurrentIndexMap(neighbor,workset,currentIndexMap)
                        if rCount and not firstTime:
                            Refine(workset.copy())
        if not updated:
            for item in sameIndexAtoms:
                workset.remove(item)

def RefineWorklist(worklist):
    workset=set(worklist)
    for item in worklist:
        if item.isComplete:
            workset.remove(item)
        else:
            item.updateNeighborIndexs()          
    Refine(workset)    

def canonize_molecule(molecule,no_isomerism=False,no_H=False,stereo=False):
    ringInfo = molecule.GetRingInfo()
    mol_conf = molecule.GetConformer()
    worklist=[atom(a,index,
                   addRDKitStereo=stereo,
                   minAtomRingSize=ringInfo.MinAtomRingSize(index),
                   coord=mol_conf.GetAtomPosition(index))
         for index,a in enumerate(molecule.GetAtoms())]
    CompleteNeighbor(worklist)
    for item in worklist:
        item.completeNeighborElements()
    for item in worklist:
        item.assignIdCode(no_H=True)
    # Initialize values
    worklist.sort(key=lambda i:i.idcode)
    workset=set(worklist)
    lastID=0
    i=0
    indexCollection=[]
    for index,item in enumerate(worklist):
        if item.idcode!=lastID:
            lastID=item.idcode
            i=index
        item.currentIndex=i
        indexCollection.append(i)
    for item in worklist:
        item.updateNeighborIndexs()
        if indexCollection.count(item.currentIndex)==1:
            item.isComplete=True
            workset.remove(item)
    # Iterative refinement
    Refine(workset)
    if not no_isomerism:
        for item in worklist:
            item.decideDoubleBondConnectedAtom()
        
        StereoAssigner(set(worklist))
        RefineWorklist(worklist)
        StereoAssigner(set(worklist),False)
        RefineWorklist(worklist)
    return worklist





