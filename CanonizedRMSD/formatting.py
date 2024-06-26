from io import StringIO
from typing import Tuple

from rdkit import Chem
from scipy import optimize
import numpy as np
import os

from CanonizedRMSD.data import *


def SequenceExchanger(f1,f2,dictionary):
    sequence=[i["canonized"] for i in sorted(dictionary,key=lambda p:p["original"])]
    substitution=[i["original"] for i in sorted(dictionary,key=lambda p:p["canonized"])]
    if dictionary[0].__contains__("item"):
        stereo=[i["item"].stereochemistry for i in sorted(dictionary,key=lambda p:p["canonized"])]
        original_stereo=[i["item"].stereochemistry for i in sorted(dictionary,key=lambda p:p["original"])]
    else:
        stereo=[0 for _ in dictionary]
        original_stereo=[0 for _ in dictionary]
    # bondstereo=[bond.GetStereo() for bond in f1.GetBonds()]
    s_file=StringIO()
    s_file.write(Chem.MolToMolBlock(f1))
    s_file.seek(0)
    content=[]
    content.append("\n     RDKit          \n\n")
    
    s_file.readline()
    s_file.readline()
    s_file.readline()
    l=s_file.readline()
    
    nAtoms=int(l[:3])
    nBonds=int(l[3:6])
    content.append(l[:27]+"  0999 V2000\n")
    Atoms=[]
    for i in range(0,nAtoms):
        Atoms.append(s_file.readline())
    
    for i in range(0,nAtoms):
        line=Atoms[substitution[i]]
        if line[30:32].strip()=="C":
            if stereo[i]==S_TYPE:
                content.append(line[:41]+"S"+line[42:])
            elif stereo[i]==R_TYPE:
                content.append(line[:41]+"R"+line[42:])
            elif stereo[i]==s_TYPE:
                content.append(line[:41]+"s"+line[42:])
            elif stereo[i]==r_TYPE:
                content.append(line[:41]+"r"+line[42:])
            else:
                content.append(line[:41]+"0"+line[42:])
        else:
            content.append(line[:41]+"0"+line[42:])

    for _ in range(0,nBonds):
        l=s_file.readline()
        atomANum=int(l[:3])
        atomBNum=int(l[3:6])
        bondOrder=int(l[6:9])
        if bondOrder!=2:
            bondStereoTag='0'
        else:
            if original_stereo[atomANum-1]==Z_TYPE and original_stereo[atomBNum-1]==Z_TYPE:
                bondStereoTag='Z'
            elif original_stereo[atomANum-1]==E_TYPE and original_stereo[atomBNum-1]==E_TYPE:
                bondStereoTag='E'
            else:
                bondStereoTag='0'
        content.append("%3s%3s"%(sequence[int(l[:3])-1]+1,sequence[int(l[3:6])-1]+1)+l[6:11]+bondStereoTag+l[12:])
    content.append("M  END\n$$$$\n")
    string=''.join(content)    # content=[]    
    if f2!=0:
        t_file=open(f2,'w')
        t_file.write(string)
        t_file.close()
    s_file.close()
    return string



