from rdkit import Chem

class atom:
    def __init__(self,a=0,OriginalIndex=0,CopyFrom=0):
        if CopyFrom==0:
            self.degree=a.GetDegree()
            self.atomicNumber=a.GetAtomicNum()
            self.attatchedHs=a.GetTotalNumHs()
            self.charge=a.GetNumRadicalElectrons()
            self.stereochemistry=int(a.GetChiralTag())  #0-Unspecified,1-,2-
            self.bonds=a.GetBonds()
            self.source=a
            self.neighbors=()
            self.neighborIndexs=[]
            self.originalIndex=OriginalIndex
            self.currentIndex=0
            self.idcode=self.__getIdCode()
            self.isComplete=False
        else:
            self.source=CopyFrom.source
            self.neighbors=()
            self.neighborIndexs=[]
            self.originalIndex=CopyFrom.originalIndex
            self.currentIndex=CopyFrom.currentIndex
            self.isComplete=CopyFrom.isComplete

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
        return idCode
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

def Refine(workset):
    worksetBackup=workset.copy()
    dictionaryForIndexs={}
    for item in workset:
        AddItem(dictionaryForIndexs,item.currentIndex,item)
    while len(workset)>0:
        # Select the partition with most instances as current partition
        currentIndexs=[item.currentIndex for item in workset]
        indexNumbers=[]
        for number in set(currentIndexs):
            indexNumbers.append({"index":number,"count":currentIndexs.count(number)})
        selectedItem=sorted(indexNumbers,key=lambda x:-x["count"])[0]
        currentIndex=selectedItem["index"]
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
            # Update neibor indexs of changed atoms' neibors and recheck them
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

def Canonizer(molecule,Branching=False):
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

    if Branching:
        collection=[]
        BranchingTieBreaking(worklist,collection)
        return collection
    else:
        OrdinaryTieBreaking(worklist)
        return worklist

def CanonizedSequenceRetriever(sourceFile,serial=False):
    mol=Chem.MolFromMolBlock(sourceFile)
    if mol==None:
        import sys
        sys.exit()
    collection=Canonizer(mol,serial)
    if serial:
        result=[]
        for canonized in collection:
            dictionary=[{"original":item.originalIndex,"canonized":item.currentIndex} for item in canonized]
            result.append(dictionary)
        return result
    else:
        dictionary=[{"original":item.originalIndex,"canonized":item.currentIndex,"item":item} for item in collection]
        return dictionary

def exitSys(start_time):
    print("Two input molecules are not identical!")
    import sys
    from time import clock
    end_time=clock()
    print("time used:%f"%(end_time-start_time))
    sys.exit(0)

def JudgeIdentity(serialA,serialB,start_time):
    for item in serialA:
        atomA=item["item"]
        atomBsource=[itemB["item"] for itemB in serialB if itemB["canonized"]==item["canonized"]]
        if len(atomBsource)==0:
            exitSys(start_time)
        atomB=atomBsource[0]
        if atomA.neighborIndexs !=atomB.neighborIndexs:
            exitSys(start_time)
    print("Two input molecules are identical!")


if __name__=='__main__':
    
    content=Chem.MolFromMolFile('testsets/test.sdf')
    collections=Canonizer(content)
    for item in collections:
        print(item.originalIndex,item.currentIndex)
    # for canonized in collections:
    #     for item in canonized:
    #         print(item.originalIndex,item.currentIndex)
        # print("===================================")


