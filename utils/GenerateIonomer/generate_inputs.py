#!/user/bin/python
import numpy as np
import itertools
import generate_system as gs
import os
import time
import utilities as utl

## function to get all unqiue combinations of two lists
def combineLists(list1,list2):
    list1 = utl.getString(list1)
    list2 = utl.getString(list2)
    if type(list1) == list:  #Convert lists into numpy arrays
        list1= np.array(list1).reshape(len(list1),1)
    elif type(list1) == str:
        list1 = np.array(list1).reshape(1,1)
    if type(list2) == list: #Convert lists into numpy arrays
        list2= np.array(list2).reshape(len(list2),1)
    elif type(list2) == str:
        list2 = np.array(list2).reshape(1,1)
    dim1 = list1.shape
    dim2 = list2.shape
    #Initialize storage
    #all_comb = np.empty([dim1[0] * dim2[0],dim1[1]+dim2[1]],dtype = object)
    all_comb = np.full([dim1[0]*dim2[0],dim1[1]+dim2[1]],'aaaaaaaaaaaaaaaa')
    #storage counter for rows
    storCounter1 = 0
    for i in range(dim1[0]):
        for j in range(dim2[0]):
            #Storage counter for columns
            storCounter2=0
            for k in range(dim1[1]):
                all_comb[storCounter1,storCounter2]= list1[i,k]
                storCounter2 +=1
            for k in range(dim2[1]):
                all_comb[storCounter1,storCounter2]= list2[j,k]
                storCounter2 +=1
            storCounter1 +=1
    # Remove any duplicate entries from duplicated inputs
    all_comb = np.unique(all_comb,axis = 0)
    return all_comb
## Function to generate the name for directory
def fileNameGen(header,propNames,propList,jobNum,listNum):
    fileName = header
    for k in range(len(propNames)):
        fileName = fileName + '_' + propertyName[k]+propList[k]
    fileName = fileName + '_' + str(jobNum+1)
    return fileName
## Function to create file with properties
def propFilesGen(propNames,propList,direcName):
    otp = open( 'overwriteProps' , 'w')
    wrt = 'direc = ' + direcName + '\n'
    otp.write(wrt)
    for k in range(len(propNames)):
        wrt = propNames[k] + ' = ' + str(propList[k]) +'\n'
        otp.write(wrt)
    return()

## Create lists of properties to vary
#Lists of functioning overwrites
  #gen_system: rho0, Na,Nb, box(change side length, cubic/square only), 'relativeSaltConc[1,2,3] (name: saltA)
    # should include but need to check: inputs that have array format
    # frac1A- sets fraction of a, keeping number of a constant

  #write_prop: spring_constant/k, chiN, chiAB, bjerrum_length
    # Special keyword- ramp- give final mulutplier for initial- ramp2-> 50 inital->100 final
    

# Notes on property naming schemes:
  # Leading zeros are used to mark values with just decimal places, so 01 is 0.1
  # p is used for strings with values before and after decimal place, so 1p5 is 1.5
     #This will overwrite leading zeros
propertyName = ['chiN','Nb','ramp']
prop1 = [20,30]
prop2 = [15,25,35]
prop3 = [1.5]


# Takes inputs for property names as lists, np arrays, single floats or strings
# Any input that produces a string longer thant 'aaaaaaaaaaaaaaaa' will fail - either make it shorter, or add more 'a's to where that squence apprears in the code (two spots)
# also might try to get it to work with lists like [0.2,0] for relativeSaltconc
# it also does not work for more than 3 properties- breaks on comb3- see if I can fix this

jobNumbers = 1    # The number of times each job is run
fileHeader = 'RepoTianrenLong' #Initial name stuff before properties in directory name

comb1 = combineLists(prop1, prop2)
comb2 = combineLists(comb1,prop3)

#np.array(utl.getString(prop1)).reshape(-1,1) #If i need to just feed it a singular list

#This needs to be manually changed to adjust
propList = comb2
## begin main loop
ln = 'Overwriting the following properties for a total of %i jobs:' % np.multiply(propList.shape[0],jobNumbers)
print(ln)
print(propertyName)
print(propList)
time.sleep(1)
for i in range(propList.shape[0]):
    for j in range(jobNumbers):
        direcName = fileNameGen(fileHeader,propertyName,propList[i,:],j,i)
        propFilesGen(propertyName,propList[i,:],direcName)
        # Call the generate system (formerly main)
        if not os.path.exists(direcName):
            os.makedirs(direcName)
        gs.genSystem()
        if os.path.exists('overwriteProps')== True:
            os.remove('overwriteProps')  #Delete file with overwrite properties
