from rdkit import Chem
import numpy as np

R_TYPE=Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW
S_TYPE=Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW
Z_TYPE=Chem.rdchem.BondStereo.STEREOZ
E_TYPE=Chem.rdchem.BondStereo.STEREOE

class atom:
    def __init__(self,a=0,OriginalIndex=0,CopyFrom=0):
        if CopyFrom==0:
            self.source=a
            self.degree=a.GetDegree()
            self.atomicNumber=a.GetAtomicNum()
            self.attatchedHs=a.GetTotalNumHs()
            self.charge=a.GetNumRadicalElectrons()
            self.stereochemistry=self.__getStereoNum()
            self.bonds=[bond.GetBondType() for bond in a.GetBonds()]
            self.neighbors=()
            self.neighborIndexs=[]
            self.originalIndex=OriginalIndex
            self.currentIndex=0
            self.idcode=self.__getIdCode()
            self.isComplete=False
            self.coordinate=molConf.GetAtomPosition(a.GetIdx())
        else:
            self.source=CopyFrom.source
            self.neighbors=()
            self.neighborIndexs=[]
            self.originalIndex=CopyFrom.originalIndex
            self.currentIndex=CopyFrom.currentIndex
            self.isComplete=CopyFrom.isComplete
            self.stereochemistry=CopyFrom.stereochemistry

    def __getIdCode(self):
        #    ***      *        *            *        *       **** 
        #atomicNum  degree  attatchedHs  charge  stereoChem  bonds(single/double/triple/aromatic)
        idCode=0
        idCode+=self.atomicNumber*1e8
        idCode+=self.degree*1e7
        idCode+=self.attatchedHs*1e6
        idCode+=self.charge*1e5
        idCode+=self.stereochemistry*1e4
        idCode+=self.bonds.count(Chem.rdchem.BondType.SINGLE)*1e3
        idCode+=self.bonds.count(Chem.rdchem.BondType.DOUBLE)*100
        idCode+=self.bonds.count(Chem.rdchem.BondType.TRIPLE)*10
        idCode+=self.bonds.count(Chem.rdchem.BondType.AROMATIC)
        return int(idCode)

    def __getStereoNum(self):
        # 0: unspecified  1: S(tetrahedron atom)  2: R(tetrahedron)  3:  E(double bond)  4:  Z(double bond)
        doubleBondList=[bond for bond in self.source.GetBonds() if bond.GetBondType()==Chem.rdchem.BondType.DOUBLE]
        if len(doubleBondList):
            doubleBond=doubleBondList[0]
            if int(doubleBond.GetStereo()):
                if doubleBond.GetStereo()==Z_TYPE:
                    return 4
                elif doubleBond.GetStereo()==E_TYPE:
                    return 3
        else:
            if self.source.HasProp("_CIPCode"):
                CIPCode=self.source.GetProp("_CIPCode")
                if CIPCode=='R':
                    return 2
                elif CIPCode=='S':
                    return 1
        return 0

    def updateNeighborIndexs(self):
        self.neighborIndexs=[item.currentIndex for item in self.neighbors]
        self.neighborIndexs.sort()

def deepcopy(worklist):
    newlist=[]
    for item in worklist:
        newlist.append(atom(CopyFrom=item))
    CompleteNeighbor(newlist)
    return newlist


def CompleteNeighbor(worklist):
    # Form neighbor network in worklist
    for item in worklist:
        item.neighbors=tuple(p for p in worklist if (p.source.GetIdx() in [q.GetIdx() for q in item.source.GetNeighbors()]))
        item.updateNeighborIndexs()

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
    selectedItem=sorted(indexNumbers,key=lambda x:-x["count"])[0]
    return selectedItem["index"]

def Refine(workset):
    dictionaryForIndexs={}
    for item in workset:
        AddItem(dictionaryForIndexs,item.currentIndex,item)
    while len(workset)>0:        
        currentIndex=SelectBestIndex(workset)
        currentPartition=dictionaryForIndexs[currentIndex].copy()
        currentPartition.sort(key=lambda i:i.neighborIndexs)
        changedIndex=[]
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
                changedIndex.append(item)
        # Decide whether there are completed atoms
        for item in currentPartition:
            if len(dictionaryForIndexs[item.currentIndex])==1 and not item.isComplete:
                item.isComplete=True
                # print(item.originalIndex,item.currentIndex,"becomes True in A") #
        # When no element in the current partition needs index change
        if len(changedIndex)==0:
            # Remove these elements from workset
            for item in currentPartition:
                workset.remove(item)
        # When changes in current partition indexs
        else:
            # Update neighbor indexs of changed atoms' neighbors and recheck them
            for item in changedIndex:
                for neighbor in item.neighbors:
                    neighbor.updateNeighborIndexs()
                    if dictionaryForIndexs.__contains__(neighbor.currentIndex):
                        for item in dictionaryForIndexs[neighbor.currentIndex]:
                            if not item.isComplete:
                                workset.add(item) 


def BranchingTieBreaking(templateWorklist,worklistCollection):
    unbrokenIndexs=[]
    for item in templateWorklist:
        if not item.isComplete:
            unbrokenIndexs.append(item.currentIndex)
    if len(unbrokenIndexs)==0:
        worklistCollection.append(templateWorklist)
    else:
        # Form dictionary for {index : Number}
        indexNumbers=[]
        for index in set(unbrokenIndexs):
            indexNumbers.append({'index':index,'number':unbrokenIndexs.count(index)})
        # Select the index with least occuring number while being as large as possible
        selectedItem=sorted(indexNumbers,key=lambda p:(p['number'],-p['index']))[0]
        selectedIndex=selectedItem['index']
        selectedNumber=selectedItem['number']
        import math
        import itertools
        permutations=list(itertools.permutations(range(selectedNumber)))
        for i in range(math.factorial(selectedNumber)):
            newList=deepcopy(templateWorklist)
            newSet=set()
            indiceDistribution=permutations[i]
            n=0
            for item in newList:
                if item.currentIndex==selectedIndex:
                    item.currentIndex+=indiceDistribution[n]
                    item.isComplete=True
                    n+=1
                elif item.currentIndex in unbrokenIndexs:
                    newSet.add(item)

            for item in newSet:
                item.updateNeighborIndexs()
            Refine(newSet)
            BranchingTieBreaking(newList,worklistCollection)

def OrdinaryTieBreaking(worklist):
    unbrokenIndexs=[]
    for item in worklist:
        if not item.isComplete:
            unbrokenIndexs.append(item.currentIndex)
    if len(unbrokenIndexs)!=0:
        # Form dictionary for {index : Number}
        indexNumbers=[]
        for index in set(unbrokenIndexs):
            indexNumbers.append({'index':index,'number':unbrokenIndexs.count(index)})
        # Select the index with least occuring number while being as large as possible
        selectedItem=sorted(indexNumbers,key=lambda p:(p['number'],-p['index']))[0]
        selectedIndex=selectedItem['index']
        selectedNumber=selectedItem['number']
        newSet=set()
        n=0
        for item in worklist:
            if item.currentIndex==selectedIndex:
                item.currentIndex+=n
                item.isComplete=True
                n+=1
            elif item.currentIndex in unbrokenIndexs:
                newSet.add(item)
        for item in newSet:
            item.updateNeighborIndexs()
        Refine(newSet)
        OrdinaryTieBreaking(worklist)

def DecideStereo(centralAtom):
    neighborAtoms=centralAtom.neighbors
    if len(neighborAtoms)==4: # connected 4 different groups
        a,b,c,d=[np.array(neighborAtom.coordinate) for neighborAtom in sorted(neighborAtoms,key=lambda x:x.currentIndex,reverse=True)]
        # assign a,b,c,d to atomic coordinates according to their current index
        o=np.array(centralAtom.coordinate)
        ab=b-a
        ac=c-a
        od=d-o
        if np.linalg.det(np.vstack([ab,ac,od]))>0: # ab x ac · od > 0 
            return R_TYPE
        else:
            return S_TYPE
    elif len(neighborAtoms)==3:
        a,b,c=[np.array(neighborAtom.coordinate) for neighborAtom in sorted(neighborAtoms,key=lambda x:x.currentIndex,reverse=True)]
        o=np.array(centralAtom.coordinate)
        ab=b-a
        ac=c-a
        oa=a-o
        if np.linalg.det(np.vstack([ab,ac,oa]))>0:
            return S_TYPE
        else:
            return R_TYPE


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

def GetAtomConnectedByDoubleBond(atom):
    doubleBond=[bond for bond in atom.source.GetBonds() if bond.GetBondType()==Chem.rdchem.BondType.DOUBLE][0]
    doubleBondAtom=[item for item in atom.neighbors if item in (doubleBond.GetBeginAtom(),doubleBond.GetEndAtom())][0]

def StereoAssigner(workset):
    while len(workset):
        currentIndex=SelectBestIndex(workset)
        for item in workset:
            if item.currentIndex==currentIndex:
                currentAtom=item
                break
        neighborCount=len(currentAtom.neighborIndexs)
        sameIndexAtoms=[item for item in workset if item.currentIndex==currentIndex]
        updated=False
        if currentAtom.atomicNumber==6: # only C atoms need to be assigned stereo tag
            if str(currentAtom.idcode)[-3]=='1' and currentAtom.geometric==0:  # one double bond, not assigned geometric tag
                if neighborCount>1 and len(set(currentAtom.neighborIndexs))==neighborCount:
                    doubleBondAtom=GetAtomConnectedByDoubleBond(currentAtom)
                    if len(doubleBondAtom.neighborIndexs)>1 and len(set(doubleBondAtom.neighborIndexs)==len(doubleBondAtom)):
                        # both sides of the double bond have different groups
                        zCount=0
                        doubleBondAtoms=[]
                        for item in sameIndexAtoms:
                            doubleBondAtom=GetAtomConnectedByDoubleBond(currentAtom)
                            doubleBondAtoms.append(doubleBondAtom)
                            geometricTag=DecideDoubleBondGeometric(item,doubleBondAtom)
                            if geometricTag==Z_TYPE:
                                zCount+=1
                            item.geometric=geometricTag
                            doubleBondAtom.geometric=geometricTag
                        updated=True
                        for atom1,atom2 in zip(sameIndexAtoms,doubleBondAtoms):
                            workset.remove(atom1)
                            workset.remove(atom2)
                            
                            if atom1.geometric==Z_TYPE:
                                atom1.currentIndex+=len(sameIndexAtoms)-zCount
                                atom2.currentIndex+=len(sameIndexAtoms)-zCount
                                updateNeighborList=atom1.neighbors
                                updateNeighborList.append(atom2.neighbors)
                                for neighbor in updateNeighborList:
                                    neighbor.updateNeighborIndexs()
                                    if neighbor.isComplete or neighbor.stereochemistry==0:
                                        workset.add(neighbor)
                        if zCount:
                            Refine(workset)
            elif str(currentAtom.idcode)[-3]=='0' and currentAtom.stereochemistry==0:  # no double bond, not assigned stereo tag                    
                rCount=0        
                if neighborCount>2 and len(set(currentAtom.neighborIndexs))==neighborCount:
                # there is no same thing in current atom's neighbor list (implicit hydrogen included)                 
                    for item in sameIndexAtoms:
                        stereoTag=DecideStereo(item)
                        if stereoTag==R_TYPE:
                            rCount+=1
                        item.stereochemistry=stereoTag
                        item.source.SetChiralTag(stereoTag)
                    updated=True
                    for item in sameIndexAtoms:
                        workset.remove(item)
                        if item.stereochemistry==R_TYPE:
                            item.currentIndex+=len(sameIndexAtoms)-rCount  # which is sCount
                            for neighbor in item.neighbors:
                                neighbor.updateNeighborIndexs()
                                if not neighbor.isComplete or neighbor.stereochemistry==0:
                                    workset.add(neighbor)
                    if rCount:
                        Refine(workset)
        if not updated:
            for item in sameIndexAtoms:
                workset.remove(item)
        StereoAssigner(workset)
                
            

def Canonizer(molecule,Branching=False,no_isomerism=False):
    worklist=[atom(a,index) for index,a in enumerate(molecule.GetAtoms())]
    CompleteNeighbor(worklist)
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
    # if not no_isomerism:
    #     workset=set(worklist)
    #     StereoAssigner(workset)
    if Branching:
        collection=[]
        BranchingTieBreaking(worklist,collection)
        return collection
    else:
        OrdinaryTieBreaking(worklist)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
        return worklist

def CanonizedSequenceRetriever(mol,serial=False,no_isomerism=False):
    global molConf
    molConf=mol.GetConformer()
    if mol==None:
        import sys
        sys.exit()
    if not no_isomerism:
        Chem.rdmolops.AssignAtomChiralTagsFromStructure(mol)
        Chem.rdmolops.DetectBondStereoChemistry(mol,molConf)
        Chem.rdmolops.AssignStereochemistry(mol)
    collection=Canonizer(mol,serial,no_isomerism)
    if serial:
        result=[]
        for canonized in collection:
            dictionary=[{"original":item.originalIndex,"canonized":item.currentIndex,"item":item} for item in canonized]
            result.append(dictionary)
        return result
    else:
        dictionary=[{"original":item.originalIndex,"canonized":item.currentIndex,"item":item} for item in collection]
        return dictionary


def JudgeIdentity(serialA,serialB):
    for item in serialA:
        atomA=item["item"]
        atomBsource=[itemB["item"] for itemB in serialB if itemB["canonized"]==item["canonized"]]
        if len(atomBsource)==0:
            return False
        atomB=atomBsource[0]
        if not (atomA.neighborIndexs == atomB.neighborIndexs \
            and atomA.stereochemistry==atomB.stereochemistry):
            return False
    return True


if __name__=='__main__':
    
    content=Chem.MolFromMolFile('testsets/test4/1_rdkit.mol')
    CanonizedSequenceRetriever(content)


