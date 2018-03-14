import random
import sys
import formatting

def GetNumAtoms(file):
    s_file=open(file)
    s_file.readline()
    s_file.readline()
    s_file.readline()
    l=s_file.readline()
    return int(l[:3])

def CreateRandomList(number):
    sampleList=list(range(number))
    random.shuffle(sampleList)
    result=[]
    for index,num in enumerate(sampleList):
        result.append({'original':index,'canonized':num})
    return result

if len(sys.argv)==3:
    a=sys.argv[1]
    b=sys.argv[2]
    numAtoms=GetNumAtoms(a)
    formatting.SequenceExchanger(a,b,CreateRandomList(numAtoms))
    print("New file created!")