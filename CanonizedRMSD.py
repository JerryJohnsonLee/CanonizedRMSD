#!/usr/bin/env python
import sys
import argparse
import time

import numpy as np

try:
    from rdkit import Chem
    import formatting
    import main
except:
    print("RDKit module not found in your Python! Please make sure RDKit is correctly installed.")
    print("Please go to http://www.rdkit.org/docs/Install.html for more details about installation.")
    sys.exit()



def Min(myList):
    indexes=list(range(len(myList)))
    package=zip(indexes,myList)
    minimum=min(package,key=lambda x:x[1])
    if minimum[1]<1e-10:
        minimum=(minimum[0],0)
    return minimum


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

def Calculate(source1,source2,saveMediates=False,outputInterrelationship=False,no_isomerism=False,no_alignment=False,qcp=False,removeHs=False,tiebreaking=True):  
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
    molA,A=formatting.Read(source1,appending[0],file1State,removeHs)
    molB,B=formatting.Read(source2,appending[1],file2State,removeHs)
    molRA=Chem.RemoveHs(molA)
    molRB=Chem.RemoveHs(molB)
    #start_time=clock()
    contentA,_=main.CanonizedSequenceRetriever(molA,False,no_isomerism)
    contentB,unbrokenB=main.CanonizedSequenceRetriever(molB,False,no_isomerism)
    canonizedA=formatting.SequenceExchanger(molA,appending[2],contentA)
    contentRA,_=main.CanonizedSequenceRetriever(molRA,False,no_isomerism)
    contentRB,_=main.CanonizedSequenceRetriever(molRB,False,no_isomerism)
    if not main.JudgeIdentity(contentRA,contentRB):
        if saveMediates:
            formatting.SequenceExchanger(molB,appending[3],contentB)
        print("Two input molecules are not identical!")
        #print(clock()-start_time)
        sys.exit()
    print("Two input molecules are identical!")
    (ma,ea)=formatting.FormMat(canonizedA)
    if tiebreaking:
        minRmsd,canonizedMinB,contentMinB=main.CanonizedSequenceRetriever(molB,True,no_isomerism,unbrokenB,ma,ea,no_alignment,qcp) 
    else:
        contentMinB=contentB
        canonizedMinB=formatting.SequenceExchanger(molB,appending[3],contentMinB)
        (mb,eb)=formatting.FormMat(canonizedMinB)
        if formatting.CheckElements(ea,eb):
            if (no_alignment):
                minRmsd=np.linalg.norm(ma-mb)/np.sqrt(ma.shape[0])
            elif (qcp==False):
                minRmsd=formatting.kabsch_rmsd(ma,mb,no_alignment=no_alignment)
            else:
                MA=ma.A
                MB=mb.A
                minRmsd=formatting.qcp_rmsd(MA,MB)
    print('RMSD='+str(minRmsd))
    if saveMediates:
        canonizedB=formatting.SequenceExchanger(molB,appending[3],contentMinB)
        _,transition,rotation,coords=formatting.kabsch_rmsd(formatting.FormMat(canonizedA)[0] \
        ,formatting.FormMat(canonizedB)[0],True)
        with open("conversion_matrices.log",'w') as f:
            s="Transition Matrix:\n"+str(transition) \
        +"\nRotation Matrix:\n"+str(rotation)+"\nTransformed Coordinates for molecule 2:\n"+str(coords)
            f.write(s)
    if outputInterrelationship:
        OutputInterrelationship(GetInterrelationship(contentA,contentMinB),A,B)

def GetConversion(molA,molB):
    contentA,_=main.CanonizedSequenceRetriever(molA,no_isomerism=True)
    contentB,unbrokenB=main.CanonizedSequenceRetriever(molB,no_isomerism=True)
    canonizedA=formatting.SequenceExchanger(molA,0,contentA)
    if not main.JudgeIdentity(contentA,contentB):
        return None
    contentBseries=main.CanonizedSequenceRetriever(molB,True,unbrokenMolecule=unbrokenB)
    rmsdCollection=[]
    for contentB in contentBseries:        
        canonizedB=formatting.SequenceExchanger(molB,0,contentB)
        ma,_=formatting.FormMat(canonizedA)
        mb,_=formatting.FormMat(canonizedB)
        rmsdCollection.append(formatting.kabsch_rmsd(ma,mb))
    
    if len(rmsdCollection)!=0:
        
        minIndex,_=Min(rmsdCollection)
        canonizedB=formatting.SequenceExchanger(molB,0,contentBseries[minIndex])
        rmsd,transition,rotation,_=formatting.kabsch_rmsd(formatting.FormMat(canonizedA)[0] \
            ,formatting.FormMat(canonizedB)[0],True)
        sortedA=[i["canonized"] for i in sorted(contentA,key=lambda p:p["original"])]
        sortedB=[i["original"] for i in sorted(contentBseries[minIndex],key=lambda p:p["canonized"])]
        mapping=dict()
        for indexA,canonized in enumerate(sortedA):
            mapping[indexA]=sortedB[canonized]  # {uncanonizedA: uncanonizedB}
        return mapping,rmsd,transition,rotation



if __name__=="__main__":
    global file1State
    global file2State

    parser=argparse.ArgumentParser( \
    description="to calculate the RMSD of two molecules after canonizing them. \
    \n\nsupported file types:\n   .mol | .sdf | .rxn | .mol2 | .ml2  | .pdb\n", \
    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("file1")
    parser.add_argument("file2")
    parser.add_argument("-s","--save",action="store_true",help="save intermediate results")
    parser.add_argument("-m","--mapping",action="store_true",help="output atom mapping relationship with two molecules")
    parser.add_argument('-i',"--ignore_isomerism",action="store_true",help="ignore geometric and stereometric isomerism when canonizing")
    parser.add_argument('-na',"--no_alignment",action="store_true",help="do not apply molecule alignment Kabsch or QCP algorithm when calculating RMSD")
    parser.add_argument('-alg', "--algorithm", choices=["Kabsch","QCP"],default="Kabsch",help="algorithm to calculate RMSD, default to Kabsch")
    parser.add_argument('-r',"--removeHs",action="store_true",help="remove H atoms")
    parser.add_argument('-at',"--arbitrary_tiebreaking",action="store_true",help="apply an arbitrary tiebreaking, namely not performing branching tiebreaking")

    args=parser.parse_args()  

    # check file format and validity
    file1State = main.CheckValidity(args.file1)
    file2State = main.CheckValidity(args.file2)
    if file1State >= 0 and file2State >= 0:
        if file1State == 0 or file2State == 0:
            print("\nWarning: Unknown file type! \n\nType -h to get supported file types.\n")
    else:
        print("\nError: Unsupported file type! \n\nType -h to get supported file types.\n")
        sys.exit(0)

    if args.save and args.no_alignment:
        print("saving mode unavailable when alignment is not done.")
        sys.exit()
    use_qcp = (args.algorithm == "QCP")
    branching_tiebreaking = not args.arbitrary_tiebreaking
    Calculate(args.file1,args.file2,args.save,args.mapping,args.ignore_isomerism,args.no_alignment,use_qcp,args.removeHs,branching_tiebreaking)

