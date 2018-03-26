from io import StringIO

from rdkit import Chem
from numpy import *

Z_TYPE=Chem.rdchem.BondStereo.STEREOZ
E_TYPE=Chem.rdchem.BondStereo.STEREOE

def ConvertFromGaussianToRdkit(f1,f2): 
    try:  
        s_file=open(f1,'r')
    except:
        print("Cannot open file %s, please check!"%f1)
        import sys
        sys.exit()
    content=[]
    content.append("\n     RDKit          \n\n")
    s_file.readline()
    s_file.readline()
    s_file.readline()
    l=s_file.readline()

    if(len(l)==40):    # already Rdkit type
        s_file.seek(0)
        
        return s_file.read()
    content.append(l[:27]+"  0999 V2000\n")

    content.append(s_file.read())
    content.append("M  END\n$$$$\n")
    
    string=''.join(content)

    if f2!=0:             # if f2=0, only get converted content, no writing
        t_file=open(f2,'w')
        t_file.write(string)
        t_file.close()
    s_file.close()
    return string

def removeHs(content):
    s_file=StringIO()
    s_file.write(content)
    s_file.seek(0)

    content=[]
    content.append("\n     RDKit          \n\n")

    s_file.readline()
    s_file.readline()
    s_file.readline()
    l=s_file.readline()
    
    nAtoms=int(l[:3])
    Atoms=[]
    nBonds=int(l[3:6])
    sequence=[]
    index=1
    for i in range(0,nAtoms):
        line=s_file.readline()
        Atoms.append(line)
        if 'H  ' != line[31:34]:
            content.append(line)
            sequence.append(index)
            index+=1
        else:
            sequence.append(0)
    countBonds=0
    for j in range(0,nBonds):
        l=s_file.readline()
        n1=int(l[:3])
        n2=int(l[3:6])
        if sequence[n1-1]*sequence[n2-1]!=0:
            content.append("%3s%3s"%(sequence[n1-1],sequence[n2-1])+l[6:12]+'\n')
            countBonds+=1
    countAtoms=nAtoms-sequence.count(0)
    content.insert(1,"%3d%3d"%(countAtoms,countBonds)+"  0  0  0  0  0  0  0  0999 V2000\n")
    content.append("M  END\n$$$$\n")
    string=''.join(content)
    return (string,sequence)


#RMSD_module

def FormMat(address):
    if ('\n' in address):
        f=StringIO()
        f.write(address)
        f.seek(0)
        
    else:
        f=open(address)
    
    lines=[]
    points=[]
    elements=[]
    for line in f:
        lines.append(line)
    f.close()
    nAtom=int(lines[3][:3])
    for n in range(4,4+nAtom):
        l=lines[n]
        points.append([float(l[:10]),float(l[10:20]),float(l[20:30])])
        elements.append(l[31:34])
    return (matrix(points),elements)

def GetIndexList(file):
    # get the non-H atom indices and form a list
    with open(file) as f:
        line=f.readline().strip()
        while line[-4:]!="ATOM":
            line=f.readline().strip()
        index=1
        sequence=[]
        while line[-4:]!="BOND":
            line=f.readline().strip()
            if line[-1]!="H":
                sequence.append(index)
            index+=1
    del sequence[-1]
    return sequence

def CheckElements(e1,e2):
    n=len(e1)
    if (len(e2)!=n):
        print('Elements don\'t match, the input molecules are not identical!')
        return 0
    else:
        for i in range(0,n):
            if(e1[i]!=e2[i]):
                print('Elements don\'t match, the input molecules are not identical!')
                return 0
        return 1

def standardSVD(matrix):
    u,s,v=linalg.svd(matrix)
    S=zeros((u.shape[1],v.shape[0]))
    S[:len(s),:len(s)]=diag(s)
    return u,S,v

def RMSD(m1,m2,mediates=0):
    m1Center=m1.sum(axis=0)/m1.shape[0]
    m2Center=m2.sum(axis=0)/m2.shape[0]
    p=m1-m1Center
    q=m2-m2Center
    a=dot(p.T,q)
    v,s,wt=standardSVD(a)
    d=sign(linalg.det(dot(v,wt).T))
    e=eye(3)
    e[2,2]=d
    u=dot(dot(wt.T,e),v.T)     # rotation matrix
    q=dot(q,u)      # transformed q coordinates
    if mediates:
        with open(mediates,'w') as f:
            s="Transition Matrix:\n"+str(m1Center-m2Center) \
            +"\nRotation Matrix:\n"+str(u)+"\nTransformed Coordinates:\n"+str(q)
            f.write(s)
    rmsd=linalg.norm(p-q)/sqrt(m1.shape[0])
    return rmsd

#  Sequence changing module

def SequenceExchanger(f1,f2,dictionary):
    sequence=[i["canonized"] for i in sorted(dictionary,key=lambda p:p["original"])]
    substitution=[i["original"] for i in sorted(dictionary,key=lambda p:p["canonized"])]
    stereo=[i["item"].stereochemistry for i in sorted(dictionary,key=lambda p:p["canonized"])]
    bondstereo=[bond.GetStereo() for bond in f1.GetBonds()]
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
        if stereo[i]==1:
            content.append(line[:41]+"S"+line[42:])
        elif stereo[i]==2:
            content.append(line[:41]+"R"+line[42:])
        else:
            content.append(line[:41]+"0"+line[42:])

    for j in range(0,nBonds):
        l=s_file.readline()
        if bondstereo[j]==Z_TYPE:
            bondStereoTag='Z'
        elif bondstereo[j]==E_TYPE:
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

def ReadFromMol(file,appending):
    A,removedHA=removeHs(ConvertFromGaussianToRdkit(file,appending))
    molA=Chem.MolFromMolBlock(A)
    return molA,removedHA

def ReadFromMol2(file):
    molA=Chem.MolFromMol2File(source1)
    removedHA=GetIndexList(source1)
    return molA,removedHA

def Read(file,appending,fileState):
    if fileState==1:
        mol,removedH=ReadFromMol(file,appending)
    elif fileState==2:
        mol,removedH=ReadFromMol2(file)
    elif fileState==0:
        try:
            mol,removedH=ReadFromMol(file,appending)
        except:
            try:
                mol,removedH=ReadFromMol2(file)
            except:
                print("Unsupported file: %s"%file)
                import sys
                sys.exit()
    return mol,removedH

if __name__=="__main__":
    print(GetIndexList("testsets/test4/r_rdkit.mol2"))