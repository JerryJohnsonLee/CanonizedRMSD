#!/home/jerry/anaconda3/bin/python 
import sys
import argparse
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
    # state -1: cannot open  0: unrecognized file type   1: mol type   2: mol2 type
    if len(f1.split('.')) > 1:
        state1=-1
        extension=f1.split('.')[1].strip()
        if extension in ['sdf','mol','rxn']:
            state1=1
        elif extension in ['mol2','ml2']:
            state1=2
    else:
        state1=0
    if len(f2.split('.')) > 1:
        state2=-1
        extension=f2.split('.')[1].strip()
        if extension in ['sdf','mol','rxn']:
            state2=1
        elif extension in ['mol2','ml2']:
            state2=2
    else:
        state2=0
    if state1>=0 and state2>=0:
        if state1==0 or state2==0:
            print("\nWarning: Unknown file type! \n\nType -h to get supported file types.\n")
        return state1,state2
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

def Calculate(source1,source2,saveMediates=False,outputInterrelationship=False,no_isomerism=False):  
    if(saveMediates):
        if(file1State):
            address1=source1.split('.')[0]
        else:
            address1=source1
        if(file2State):
            address2=source2.split('.')[0]
        else:
            address2=source2
        appending=[address1+"_rdkit.mol",address2+"_rdkit.mol",address1+"_canonized.mol",address2+"_canonized.mol"]
    else:
        appending=[0,0,0,0]
    start_time=clock()
    molA,removedHA=Read(source1,appending[0],file1State)
    contentA=main.CanonizedSequenceRetriever(molA,False,no_isomerism)
    molB,removedHB=Read(source2,appending[1],file2State)    
    contentB=main.CanonizedSequenceRetriever(molB,False,no_isomerism)
    canonizedA=SequenceExchanger(molA,appending[2],contentA)
    if not main.JudgeIdentity(contentA,contentB):
        if saveMediates:
            SequenceExchanger(molB,appending[3],contentB)
        print("Two input molecules are not identical!")
        sys.exit()
    print("Two input molecules are identical!")
    end_time_1=clock()
    contentBseries=main.CanonizedSequenceRetriever(molB,True,no_isomerism)    
    rmsdCollection=[]
    for contentB in contentBseries:        
        canonizedB=SequenceExchanger(molB,0,contentB)
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
            canonizedB=SequenceExchanger(molB,appending[3],contentBseries[minIndex])
            RMSD(FormMat(canonizedA)[0],FormMat(canonizedB)[0],mediates="conversion_matrices")
        if outputInterrelationship:
            OutputInterrelationship(GetInterrelationship(contentA,contentBseries[minIndex]),removedHA,removedHB)
    end_time_2=clock()
    print('judging time:%f,calculation time:%f'%(end_time_1-start_time,end_time_2-start_time))

if __name__=="__main__":
    parser=argparse.ArgumentParser( \
    description="to calculate the RMSD of two molecules after canonizing them. \
    \n\nsupported file types:\n   .mol | .sdf | .rxn | .mol2 | .ml2 \n", \
    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("file1")
    parser.add_argument("file2")
    parser.add_argument("-s","--save",action="store_true",help="save intermediate results")
    parser.add_argument("-m","--mapping",action="store_true",help="output atom mapping relationship with two molecules")
    parser.add_argument('-i',"--ignore_isomerism",action="store_true",help="ignore geometric and stereometric isomerism when canonizing")

    args=parser.parse_args()  
    global file1State
    global file2State
    file1State,file2State=CheckValidity(args.file1,args.file2)
    Calculate(args.file1,args.file2,args.save,args.mapping,args.ignore_isomerism)
    #file1State=0
    #file2State=0
    #Calculate('testsets/test2/1.sdf','testsets/test2/2.sdf',saveMediates=True,no_isomerism=True)
