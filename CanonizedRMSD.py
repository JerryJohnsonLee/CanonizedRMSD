import sys
from time import clock

from rdkit import Chem
from collections import defaultdict
import operator
import sys

from formatting import *
import main

def help():
    print("Incorrect usage! Type -h to get help")
    sys.exit(0)

def Min(myList):
    indexes=list(range(len(myList)))
    package=zip(indexes,myList)
    minimum=min(package,key=lambda x:x[1])
    if minimum[1]<1e-10:
        minimum=(minimum[0],0)
    return minimum



def CheckValidity(f1,f2):
    if len(f1.split('.')) > 1:
        state1=0
        if f1.split('.')[1].strip() in ['sdf','mol','rxn']:
            state1=1
    else:
        state1=2
    if len(f2.split('.')) > 1:
        state2=0
        if f2.split('.')[1].strip() in ['sdf','mol','rxn']:
            state2=1
    else:
        state2=2
    if state1!=0 and state2!=0:
        if state1==2 or state2==2:
            print("\nWarning: Unknown file type! \n\nType -h to get supported file types.\n")
        return 1
    else:
        print("\nError: Unsupported file type! \n\nType -h to get supported file types.\n")
        sys.exit(0)

def GetInterrelationship(canonizedCollection1,canonizedCollection2):
    result=[]
    for item in canonizedCollection1:
        canonizedIndex=item["canonized"]

        index2=[atom["original"] for atom in canonizedCollection2 if atom["canonized"]==canonizedIndex][0]
        result.append((item["original"],index2))
    return result

def OutputInterrelationship(collection,sequenceA,sequenceB):
    print("="*40)
    print("Index in File 1 || Index in File 2")
    for item in collection:
        print(str(sequenceA.index(item[0]+1)+1).center(16)+"  "+str(sequenceB.index(item[1]+1)+1).center(16))

def Calculate(source1,source2,saveMediates=False,outputInterrelationship=False):  
    if(saveMediates):
        if(len(source1.split('.'))==2):
            address1=source1.split('.')
        elif(len(source1.split('.'))==1):
            address1=[source1,'']
        if(len(source2.split('.'))==2):
            address2=source2.split('.')
        elif(len(source2.split('.'))==1):
            address2=[source2,'']
        appending=[address1[0]+"_rdkit."+address1[1],address2[0]+"_rdkit."+address2[1],address1[0]+"_canonized."+address1[1],address2[0]+"_canonized."+address2[1]]
    else:
        appending=[0,0,0,0]
    start_time=clock()
    A,removedHA=removeHs(ConvertFromGaussianToRdkit(source1,appending[0]))
    B,removedHB=removeHs(ConvertFromGaussianToRdkit(source2,appending[1]))
    contentA=main.CanonizedSequenceRetriever(A)
    contentB=main.CanonizedSequenceRetriever(B)
    main.JudgeIdentity(contentA,contentB,start_time)
    end_time_1=clock()
    contentBseries=main.CanonizedSequenceRetriever(B,True)
    canonizedA=SequenceExchanger(A,appending[2],contentA)
    rmsdCollection=[]
    for contentB in contentBseries:        
        canonizedB=SequenceExchanger(B,0,contentB)
        (ma,ea)=FormMat(canonizedA)
        (mb,eb)=FormMat(canonizedB)
        if CheckElements(ea,eb):
            rmsdCollection.append(RMSD(ma,mb))
        else:
            break
    if len(rmsdCollection)!=0:
        minIndex,minimum=Min(rmsdCollection)
        print('RMSD='+str(minimum))
        if saveMediates:
            canonizedB=SequenceExchanger(B,appending[3],contentBseries[minIndex])
            RMSD(FormMat(canonizedA)[0],FormMat(canonizedB)[0],mediates="conversion_matrices")
        if outputInterrelationship:
            OutputInterrelationship(GetInterrelationship(contentA,contentBseries[minIndex]),removedHA,removedHB)
    end_time_2=clock()
    print('judging time:%f,calculation time:%f'%(end_time_1-start_time,end_time_2-start_time))

if __name__=="__main__":
    if len(sys.argv)==3:   
        a=sys.argv[1]
        b=sys.argv[2]
        if CheckValidity(a,b):
            Calculate(a,b)
    elif len(sys.argv)==4:
        a=sys.argv[1]
        b=sys.argv[2] 
        c=sys.argv[3]
        saveMediates=False
        outputRelations=False
        if c[0]=='-':
            if 'S' in c.upper():
                saveMediates=True
            if 'O' in c.upper():
                outputRelations=True
            Calculate(a,b,saveMediates,outputRelations)
        else:
            help()
    elif len(sys.argv)>4:
        print("Too many input arugments! Type -h to get help")
        sys.exit(0)
    elif len(sys.argv)==2:
        if sys.argv[1].upper()=='-H':
            print("\nCanonizedRMSD: \n   To calculate the RMSD of two molecules after canonizing them.\n\nUsage:\n   CanonizedRMSD.py File1 File2 [Options]\n\nOptions:\n   -s    Saving Intermediates (Not saving by default)\n   -o    Output interrelationship of corresponding atoms in File 1 and File 2 (Not outputting by default)\n\nSupported file types:\n   .mol | .sdf | .rxn\n")
        else:
            help()
    else:
        help()

