from rdkit import Chem
import numpy as np
import formatting
import os
import sys
import time

R_TYPE=4
S_TYPE=3
Z_TYPE=6
E_TYPE=5
r_TYPE=2
s_TYPE=1
AROMATIC=7
DOUBLE_BOND=Chem.rdchem.BondType.DOUBLE
AROMATIC_BOND=Chem.rdchem.BondType.AROMATIC

class atom:
    def __init__(self,a=None,OriginalIndex=0,CopyFrom=None):
        if CopyFrom==None:
            self.source=a
            self.degree=a.GetDegree()
            self.atomicNumber=a.GetAtomicNum()
            self.attatchedHs=a.GetTotalNumHs(includeNeighbors=True)
            self.charge=a.GetNumRadicalElectrons()
            self.isotope=a.GetIsotope()
            self.stereochemistry=0
            self.bonds=[bond.GetBondType() for bond in a.GetBonds()]
            self.doubleBondConnectAtom=None
            self.neighbors=()
            self.neighborIndexs=[]
            self.neighborElements=()
            self.originalIndex=OriginalIndex
            self.currentIndex=0
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
            self.doubleBondConnectAtom=CopyFrom.doubleBondConnectAtom
            self.idcode=CopyFrom.idcode

    def assignIdCode(self, loose_criterion=False):
        # new criterion:
        #    ***                ***                  *       *       *        *
        # atomicNum    largestNeighborElement   connection  n_Hs  isotope   charge
        if loose_criterion:
            connection=self.degree
        else:
            connection=self.degree+self.bonds.count(Chem.rdchem.BondType.DOUBLE) \
            +self.bonds.count(Chem.rdchem.BondType.TRIPLE)*2 \
            +self.bonds.count(Chem.rdchem.BondType.AROMATIC)*0.5
        largestNeighborElement=max([item.atomicNumber for item in self.neighbors])
        idCode=0
        idCode+=self.atomicNumber*1e7
        idCode+=largestNeighborElement*1e4
        idCode+=connection*1e3
        idCode+=self.attatchedHs*100
        idCode+=self.isotope*10
        if not loose_criterion:
            idCode+=(self.charge+4)   # to avoid problem with negative charges
        self.idcode=int(idCode)

    def completeNeighborElements(self):
        self.neighborElements=tuple(sorted([neighborAtom.source.GetAtomicNum() \
     for neighborAtom in self.neighbors],reverse=True))   # a tuple is used to make it hashable

    def updateNeighborIndexs(self):
        self.neighborIndexs=[item.currentIndex for item in self.neighbors]
        self.neighborIndexs.sort(reverse=True)

    def decideDoubleBondConnectedAtom(self):
        if self.atomicNumber==6 and DOUBLE_BOND in self.bonds:
            doubleBond=[bond for bond in self.source.GetBonds() if bond.GetBondType()==DOUBLE_BOND][0]
            doubleBondAtom=[item for item in self.neighbors \
                if item.source.GetIdx() in (doubleBond.GetBeginAtomIdx(),doubleBond.GetEndAtomIdx())][0]
            self.doubleBondConnectAtom=doubleBondAtom


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
    selectedItem=sorted(indexNumbers,key=lambda x:(-x["count"],-x["index"]))[0]
    return selectedItem["index"]

def Refine(workset):
    dictionaryForIndexs={}
    for item in workset:
        AddItem(dictionaryForIndexs,item.currentIndex,item)
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
                                workset.add(item) 


def BranchingTieBreaking(molB,templateWorklist,ma,ea,no_alignment,qcp):
    unbrokenIndexs=[]
    for item in templateWorklist:
        if not item.isComplete:
            unbrokenIndexs.append(item.currentIndex)
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
    minRmsd=-1
    canonizedMinB=None
    contentMinB=None
    List=None
    for i in range(math.factorial(selectedNumber)):
        newList=deepcopy(templateWorklist)
        newSet=set()
        indiceDistribution=permutations[i]
        n=0
        for j in newList:
            if j.currentIndex==selectedIndex:
                j.currentIndex+=indiceDistribution[n]
                j.isComplete=True
                n+=1
            elif j.currentIndex in unbrokenIndexs:
                newSet.add(j)
        for j in newSet:
            j.updateNeighborIndexs()
        Refine(newSet)
        contentB=[{"original":k.originalIndex,"canonized":k.currentIndex,"item":k} for k in newList]
        canonizedB=formatting.SequenceExchanger(molB,0,contentB)
        (mb,eb)=formatting.FormMat(canonizedB)
        if formatting.CheckElements(ea,eb):
            if (no_alignment):
                rmsd=np.linalg.norm(ma-mb)/np.sqrt(ma.shape[0])
            elif (qcp==False):
                rmsd=formatting.kabsch_rmsd(ma,mb,no_alignment=no_alignment)
            else:
                MA=ma.A
                MB=mb.A
                rmsd=formatting.qcp_rmsd(MA,MB)
        else:
            sys.exit()
        if minRmsd==-1 or minRmsd>rmsd:
            minRmsd=rmsd
            canonizedMinB=canonizedB
            contentMinB=contentB
            List=newList
    complete=True
    for item in newList:
        if not item.isComplete:
            complete=False
            break
    if complete:
        return minRmsd,canonizedMinB,contentMinB
    else:
        return BranchingTieBreaking(molB,List,ma,ea,no_alignment,qcp)

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
        # selectedNumber=selectedItem['number']
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


def StereoAssigner(workset,firstTime=True):
    while len(workset):
        currentIndex=SelectBestIndex(workset)
        for item in workset:
            if item.currentIndex==currentIndex:
                currentAtom=item
                break
        neighborCount=len(currentAtom.neighborIndexs)
        sameIndexAtoms=[item for item in workset if item.currentIndex==currentIndex] # find all atoms with same index
        updated=False
        if currentAtom.stereochemistry==0: # check if the atom has not been assigned stereo tag, start assigning
            if currentAtom.bonds.count(DOUBLE_BOND)==1 :  # one double bond, not assigned geometric tag
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
                                            workset.add(neighbor)
                            if zCount:
                                Refine(workset.copy())
            elif currentAtom.bonds.count(DOUBLE_BOND)==0 and currentAtom.bonds.count(AROMATIC_BOND)==0:  # no double bond or aromatic bond                  
                rCount=0        
                if neighborCount>2 and len(set(currentAtom.neighborIndexs))==neighborCount:
                # there is no same thing in current atom's neighbor list (implicit hydrogen included)                 
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
                        for neighbor in item.neighbors:
                            neighbor.updateNeighborIndexs()
                            if not neighbor.isComplete or neighbor.stereochemistry==0:
                                if not firstTime:
                                    workset.add(neighbor)
                    if rCount and not firstTime:
                        Refine(workset.copy())
        if not updated:
            for item in sameIndexAtoms:
                workset.remove(item)

                
            

def Canonizer(molecule,no_isomerism=False):
    worklist=[atom(a,index) for index,a in enumerate(molecule.GetAtoms())]
    CompleteNeighbor(worklist)
    for item in worklist:
        item.completeNeighborElements()
    for item in worklist:
        item.assignIdCode()
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
        workset=set(worklist)
        for item in worklist:
            if item.isComplete:
                workset.remove(item)
        Refine(workset)
        StereoAssigner(set(worklist),False)
    return worklist

def CanonizedSequenceRetriever(mol,serial=False,no_isomerism=False,unbrokenMolecule=None,ma=None,ea=None,no_alignment=False,qcp=False):
    if not unbrokenMolecule:
        global molConf
        if mol==None:
            import sys
            sys.exit()
        molConf=mol.GetConformer()
        unbrokenMolecule=Canonizer(mol,no_isomerism)
    if serial:
        return BranchingTieBreaking(mol,unbrokenMolecule,ma,ea,no_alignment,qcp)
    else:
        collection=deepcopy(unbrokenMolecule)
        OrdinaryTieBreaking(collection)
        dictionary=[{"original":item.originalIndex,"canonized":item.currentIndex,"item":item} for item in collection]
        return dictionary,unbrokenMolecule


def JudgeIdentity(serialA,serialB):
    for item in serialA:
        atomA=item["item"]
        atomBsource=[itemB["item"] for itemB in serialB if itemB["canonized"]==item["canonized"]]
        if len(atomBsource)==0:
            return False
        atomB=atomBsource[0]
        if not (atomA.neighborIndexs == atomB.neighborIndexs \
            and atomA.stereochemistry==atomB.stereochemistry):
            print(atomA.stereochemistry)
            print(atomB.stereochemistry)
            return False
    return True

def CheckValidity(filename):
    # state -1: cannot open  0: unrecognized file type   1: mol type   2: mol2 type   3: pdb type
    f1 = os.path.basename(filename)
    if len(f1.split('.')) > 1:
        state = -1
        extension=f1.split('.')[1].strip().lower()
        if extension in ['sdf','mol','rxn']:
            state = 1
        elif extension in ['mol2','ml2']:
            state = 2
        elif extension in ['pdb']:
            state = 3
    else:
        state = 0
    return state


if __name__=='__main__':
    
    content=Chem.MolFromMolFile('testsets/test4/1_rdkit.mol')
    CanonizedSequenceRetriever(content)
