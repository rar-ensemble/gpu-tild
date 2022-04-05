import numpy as np
import os, sys
from numpy.random import rand, normal

def addSalt(box,Dim, saltNumber, saltIndex, saltCharge, moleculeIndex):
    #generate initial variables
    x = np.zeros((saltNumber*2,Dim),'d') 
    tp = np.zeros((saltNumber*2),'i') 
    ql = np.zeros((saltNumber*2),'i')
    mID = np.zeros((saltNumber*2),'i')
    ind = 0       #running index counter
    for i in range(2):  #Loop for +/- charge
        for j in range(saltNumber): #Loop for number of salt
            for k in range(Dim):   #x position for each dimension
                x[ind,k] = rand()*box[k]
            ql[ind] = saltCharge * np.power(-1,i)   #Charge
            tp[ind] = saltIndex + i                 #Type
            mID[ind] = moleculeIndex + ind            # Molecule index number
            ind += 1
    return x, tp, ql, mID

def getValues(inStr):
    #This converts string of property number into acual number
    # It also converts leading zeros into decimal points according to my naming convention
    # This is used to convert my file/directory names (which can't use decimal points) to the actual decimal
    # they represent. Example 01 is 0.1, 001 is 0.01. First 0 is before ., everything else is after
    # It also used p for a decimal point with higher priority
    convertOngoing = True
    # Check if there is a 'p' in the string
    if 'p' in inStr:
        inStr2 = inStr.replace('p','.')
        inStr = inStr2
        convertOngoing = False
    convertCounter = 0;
    divisor = 1;
    while convertOngoing:   #This converts leading zeros into the appropriate number of decimal points
        if inStr[convertCounter]== '0':
            divisor *=10;
            convertCounter += 1;
            if convertCounter == len(inStr):
                print('zero properties do not work')
                convertOngoing =False
        else:
            convertOngoing = False
    outFlt = float(inStr)/divisor
    return outFlt
def getString(inpt):
    # This converts an interger or a float into the string naming convention used for my file/directory names
    # will accept sigular values, lists, or np.arrays, and will return the same type
    if type(inpt) == float or type(inpt) == int:
        otp = str(inpt)
        if '.' in otp:
            otp2 = otp.replace('.','p')
            otp = otp2
    elif type(inpt) == str:
        otp = inpt
        if '.' in otp:
            otp2 = otp.replace('.','p')
            otp = otp2
    elif type(inpt) == list:
        otp = []
        for i in range(len(inpt)):
            if type(inpt[i]) == float or type(inpt[i]) == int:
                otp.append(str(inpt[i]))
            elif type(inpt[i]) == str:
                otp.append(inpt[i])
            else:
                quit('unsupported data type in inputs')
            # Replace decimal points with (p)s
            if '.' in otp[i]:
                    otp2 = otp[i].replace('.','p')
                    otp[i] = otp2
    elif type(inpt) == np.ndarray:
        dim = inpt.shape
        otp = np.full(dim,'aaaaaaaaaaaaaaaa')
        for i in range(dim[0]):
            for j in range(dim[1]):
                if type(inpt[i,j]) == float or type(inpt[i,j]) == int:
                    otp[i,j] = str(inpt)
                elif type(inpt[i,j])==str or type(inpt[i,j])== np.str_:
                    otp[i,j]=inpt[i,j]
                else:
                    quit('Unsupported data type in inputs')
                if '.' in otp[i,j]:
                    otp2 = otp[i].replace('.','p')
                    otp[i] = otp2
    else:
        quit('Unsupported Data type in inputs')
    return otp

def periodicBC1(r,box): #check PBCs for one dimension
    if ( r > box ):
        r = r - box    # check PBC
    elif (r < 0.0 ):
        r = r + box   # check PBC
    return r

def periodicBC(r,box,Dim):
    for j in range(Dim):
        r[0,j] = perodicBC1(r[0,j],box[j],Dim[j])
    return r
