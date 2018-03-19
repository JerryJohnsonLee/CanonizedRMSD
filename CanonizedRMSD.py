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

def Calculate(source1,source2,saveMediates=False,outputInterrelationship=False,no_isomerism=False):  
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
    contentA,molA=main.CanonizedSequenceRetriever(A,False,no_isomerism)
    contentB,molB=main.CanonizedSequenceRetriever(B,False,no_isomerism)
    canonizedA=SequenceExchanger(molA,appending[2],contentA)
    if not main.JudgeIdentity(contentA,contentB):
        if saveMediates:
            SequenceExchanger(molB,appending[3],contentB)
        print("Two input molecules are not identical!")
        sys.exit()
    print("Two input molecules are identical!")
    end_time_1=clock()
    contentBseries=main.CanonizedSequenceRetriever(B,True,no_isomerism)    
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
    \n\nsupported file types:\n   .mol | .sdf | .rxn\n", \
    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("file1")
    parser.add_argument("file2")
    parser.add_argument("-s","--save",action="store_true",help="save intermediate results")
    parser.add_argument("-m","--mapping",action="store_true",help="output atom mapping relationship with two molecules")
    parser.add_argument('-i',"--no_isomerism",action="store_true",help="do not consider geometric and stereometric isomerism when canonizing")

    # args=parser.parse_args()
    # if CheckValidity(args.file1,args.file2):
    #     Calculate(args.file1,args.file2,args.save,args.mapping,args.no_isomerism)

    Calculate('testsets/test4/1.mol','testsets/test4/2.mol',saveMediates=True)
   
