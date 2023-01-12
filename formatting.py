from io import StringIO
from typing import Tuple

from rdkit import Chem
from scipy import optimize
import numpy as np
import os

R_TYPE=4
S_TYPE=3
Z_TYPE=6
E_TYPE=5
r_TYPE=2
s_TYPE=1
AROMATIC=7

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

def removeHs(content,removeH):
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
    for _ in range(0,nAtoms):
        line=s_file.readline()
        Atoms.append(line)
        if 'H  ' != line[31:34] or removeH==False:
            content.append(line)
            sequence.append(index)
            index+=1
        else:
            sequence.append(0)
    countBonds=0
    for _ in range(0,nBonds):
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

def center_of_geometry(coordinates: np.ndarray) -> np.ndarray:
    #Center of geometry.
    assert coordinates.shape[1] == 3
    return np.mean(coordinates, axis=0)

def center(coordinates: np.ndarray) -> np.ndarray:
    #Center coordinates.
    return coordinates - center_of_geometry(coordinates)

def M_mtx(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    #Compute inner product between coordinate matrices.
    return B.T @ A

def K_mtx(M):
    #Compute symmetric key matrix.
    assert M.shape == (3, 3)
    S_xx = M[0, 0]
    S_xy = M[0, 1]
    S_xz = M[0, 2]
    S_yx = M[1, 0]
    S_yy = M[1, 1]
    S_yz = M[1, 2]
    S_zx = M[2, 0]
    S_zy = M[2, 1]
    S_zz = M[2, 2]
    # p = plus, m = minus
    S_xx_yy_zz_ppp = S_xx + S_yy + S_zz
    S_yz_zy_pm = S_yz - S_zy
    S_zx_xz_pm = S_zx - S_xz
    S_xy_yx_pm = S_xy - S_yx
    S_xx_yy_zz_pmm = S_xx - S_yy - S_zz
    S_xy_yx_pp = S_xy + S_yx
    S_zx_xz_pp = S_zx + S_xz
    S_xx_yy_zz_mpm = -S_xx + S_yy - S_zz
    S_yz_zy_pp = S_yz + S_zy
    S_xx_yy_zz_mmp = -S_xx - S_yy + S_zz
    return np.array(
        [
            [S_xx_yy_zz_ppp, S_yz_zy_pm, S_zx_xz_pm, S_xy_yx_pm],
            [S_yz_zy_pm, S_xx_yy_zz_pmm, S_xy_yx_pp, S_zx_xz_pp],
            [S_zx_xz_pm, S_xy_yx_pp, S_xx_yy_zz_mpm, S_yz_zy_pp],
            [S_xy_yx_pm, S_zx_xz_pp, S_yz_zy_pp, S_xx_yy_zz_mmp],
        ]
    )


def coefficients(M: np.ndarray, K: np.ndarray) -> Tuple[float, float, float]:
    #Compute quaternion polynomial coefficients.
    c2 = -2 * np.trace(M.T @ M)
    c1 = -8 * np.linalg.det(M)  # TODO: Slow?
    c0 = np.linalg.det(K)  # TODO: Slow?
    return c2, c1, c0


def lambda_max(Ga: float, Gb: float, c2: float, c1: float, c0: float) -> float:
    #Find largest root of the quaternion polynomial.
    def P(x):
        #Quaternion polynomial
        return x ** 4 + c2 * x ** 2 + c1 * x + c0
    def dP(x):
        #Fist derivative of the quaternion polynomial
        return 4 * x ** 3 + 2 * c2 * x + c1
    x0 = (Ga + Gb) * 0.5
    lmax = optimize.newton(P, x0, fprime=dP)
    return lmax


def _lambda_max_eig(K: np.ndarray) -> float:
    #Find largest eigenvalue of :math:`K`.
    e, _ = np.linalg.eig(K)
    return max(e)


def qcp_rmsd(a: np.ndarray, b: np.ndarray, atol: float = 1e-9) -> float:
    #Compute RMSD using the quaternion polynomial method.
    
    A = center(a)
    B = center(b)

    assert A.shape == B.shape

    N = A.shape[0]

    Ga = np.trace(A.T @ A)
    Gb = np.trace(B.T @ B)

    M = M_mtx(A, B)
    K = K_mtx(M)

    c2, c1, c0 = coefficients(M, K)

    try:
        # Fast calculation of the largest eigenvalue of K as root of the characteristic
        # polynomial.
        l_max = lambda_max(Ga, Gb, c2, c1, c0)
    except RuntimeError:  # Newton method fails to converge; see GitHub Issue #35
        # Fallback to (slower) explicit calculation of the largest eigenvalue of K
        l_max = _lambda_max_eig(K)

    s = Ga + Gb - 2 * l_max

    if abs(s) < atol:  # Avoid numerical errors when Ga + Gb = 2 * l_max
        rmsd = 0.0
    else:
        rmsd = np.sqrt(s / N)

    return rmsd

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
    return (np.matrix(points),elements)

def GetIndexList(file,removeHs):
    # get the non-H atom indices and form a list
    with open(file) as f:
        line=f.readline().strip()
        while line[-4:]!="ATOM":
            line=f.readline().strip()
        index=1
        sequence=[]
        while line[-4:]!="BOND":
            line=f.readline().strip()
            if (removeHs):
                if line[-1]!="H":
                    sequence.append(index)
            else:
                sequence.append(index)
            index+=1
    del sequence[-1]
    return sequence

def GetNonHydrogenIndices(molecule):
    # get the non-H atom indices and form a list
    sequence=[]
    for atom in molecule.GetAtoms():
        if atom.GetSymbol()!="H":
            sequence.append(atom.GetIdx())
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
    u,s,v=np.linalg.svd(matrix)
    S=np.zeros((u.shape[1],v.shape[0]))
    S[:len(s),:len(s)]=np.diag(s)
    return u,S,v

def kabsch_rmsd(m1,m2,output=False,no_alignment=False):
    if no_alignment:
        return np.linalg.norm(m1-m2)/np.sqrt(m1.shape[0])
    else:
        m1Center=m1.sum(axis=0)/m1.shape[0]
        m2Center=m2.sum(axis=0)/m2.shape[0]
        p=m1-m1Center
        q=m2-m2Center
        a=np.dot(p.T,q)
        v,s,wt=standardSVD(a)
        d=np.sign(np.linalg.det(np.dot(v,wt).T))
        e=np.eye(3)
        e[2,2]=d
        u=np.dot(np.dot(wt.T,e),v.T)     # rotation matrix
        q=np.dot(q,u)      # transformed q coordinates
        rmsd=np.linalg.norm(p-q)/np.sqrt(m1.shape[0])
        if output:
            return rmsd,(m2Center-m1Center),u,q+m1Center
    
        return rmsd

#  Sequence changing module

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

def ReadFromMol(file,appending,removeH):
    convertedMolBlock=ConvertFromGaussianToRdkit(file,appending)
    molA=Chem.MolFromMolBlock(convertedMolBlock,sanitize=False)
    A=GetNonHydrogenIndices(molA)
    return molA,A

def ReadFromMol2(file,removeH):
    molA=Chem.MolFromMol2File(file,sanitize=False)
    A=GetNonHydrogenIndices(molA)
    return molA,A

def ReadFromMol3(file,appending,removeH):
    molA = Chem.MolFromPDBFile(file,sanitize=False)
    A=GetNonHydrogenIndices(molA)
    return molA,A

def Read(file,appending,fileState,removeH,aromatize=False,assignRDKitStereo=False):
    if fileState==1:
        mol,M=ReadFromMol(file,appending,removeH)
    elif fileState==2:
        mol,M=ReadFromMol2(file,removeH)
    elif fileState==3:
        mol,M=ReadFromMol3(file,appending,removeH)
    elif fileState==0:
        try:
            mol,M=ReadFromMol(file,appending,removeH)
        except:
            try:
                mol,M=ReadFromMol2(file,removeH)
            except:
                print("Unsupported file: %s"%file)
                import sys
                sys.exit()
    
    # partial sanitization to make sure necessary molecular properties are calculated
    mol.UpdatePropertyCache(strict=False)
    sanitizeOpts=Chem.SanitizeFlags.SANITIZE_FINDRADICALS| \
                 Chem.SanitizeFlags.SANITIZE_KEKULIZE| \
                 Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION| \
                 Chem.SanitizeFlags.SANITIZE_SYMMRINGS| \
                 Chem.SanitizeFlags.SANITIZE_ADJUSTHS
    if aromatize:
        sanitizeOpts=sanitizeOpts|Chem.SanitizeFlags.SANITIZE_SETAROMATICITY| \
                                  Chem.SanitizeFlags.SANITIZE_SETCONJUGATION

    sanitization_result = Chem.SanitizeMol(mol,sanitizeOps=sanitizeOpts,
                                     catchErrors=True)
    if sanitization_result != Chem.SanitizeFlags.SANITIZE_NONE:
        print("Cannot sanitize file: %s"%file)
        import sys
        sys.exit()

    # assign rdkit stereochemistry tags to be used for initializing idCode (when necessary)
    if assignRDKitStereo:
        Chem.AssignStereochemistryFrom3D(mol)
    return mol,M


if __name__=="__main__":
    print(GetIndexList("testsets/test4/r_rdkit.mol2"))
