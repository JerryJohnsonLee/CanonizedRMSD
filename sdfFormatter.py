import sys
import formatting

def help():
    print("Incorrect usage! Type -h to get help")
    sys.exit(0)




def CheckValidity(f1,f2):
    if len(f1.split('.')) > 1 and f1.split('.')[1].strip() == 'sdf':
        if len(f2.split('.')) > 1 and f2.split('.')[1].strip() == 'sdf':
            return 1
        else:
            print(".sdf format file needed for the target file")
    else:
        print(".sdf format file needed for the source file")
    sys.exit(0)

if len(sys.argv)==3:
    
    a=sys.argv[1]
    b=sys.argv[2]
    if CheckValidity(a,b):
        formatting.ConvertFromGaussianToRdkit(a,b)
        

elif len(sys.argv)==4:
    a=sys.argv[1]
    b=sys.argv[2] 
    c=sys.argv[3]
    if c.upper()  =='-R':
        if CheckValidity(a,b):
            formatting.ConvertFromGaussianToRdkit(a,b)
    elif c.upper() =='-G':
        if CheckValidity(a,b):
            formatting. ConvertFromRdkitToGaussian(a,b)
    else:
        help()
elif len(sys.argv)>4:
    print("Too many input arugments! Type -h to get help")
    sys.exit(0)
elif len(sys.argv)==2:
    if sys.argv[1].upper()=='-H':
        print("sdfformatting: \n   to format a .sdf file to a certain type\n\nUsage:\n   sdfformatting.py sourceFile targetFile [direcction]\n\nDirection:\n   -r    From Gaussian to Rdkit\n   -g    From Rdkit to Gaussian")
    else:
        help()
else:
    help()

