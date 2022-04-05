#!/home/rrig/Install/anaconda3/bin/python3
import numpy as np
import os, sys
import utilities as utl
from numpy.random import rand, normal

def write_prop(propertiesName, Dim, name, rho0, box, bl , tp, isIonomer, drudeOscillator, do_charge,saltConc,N):
    # User defined properties not defined in main
    max_steps = 3000000                             # Number of steps simulation runs
    log_freq = 5000                                 # Frequency of writing to data.dat
    binary_freq = 5000                             # Frequency of writing to binary data files
    delt = 0.005                                    # Size of time step
    pmeorder = 1                                    # Order of spline used to map particle density to grid [1,2,3,4]
    integrator = 'GJF'                              # Integrator used (EM or GJF)
    integratorGroupName = 'all'                     # Group name for integrator

    readResumeLocation = 'none'            # Where to read resume file - set to 'none' if no read_resume

    nodeRatio = 81/25
    nodeCount = [125,125,125] #np.floor(nodeRatio*box)
    # make it odd
    for i1 in range(len(nodeCount)):
        if (nodeCount[i1] % 2) == 0:
            nodeCount[i1] += 1
    
    kappa = 10                                      # Kappa for interaction (what is it)
    spring_constant = 135 #90                            # Spring constant for Drude-Oscillator
    drude_type = 3 #3- make this auto adjust                                  # What monomer type is drude-oscilltor
    chiN = 20
    chiAB = chiN/np.sum(N)                                     # Flory huggins chi parameter
    chiCounterion = 0                               # Chi for the counterion
    std_dev = 1.0                                   # Standard Deviaiton for Gaussian
    bjerrum_length = 93#1443.586                  # Bjerrum Length
    # make it smaller
    smearing_size = 0.3                             # Charge smearing size
    traj_freq = 1000000                               # Steps that trajectory file is written to
    struc_freq = 100                                # Freqency of writing structure factor
    skip_steps = 100000                             # Steps it skips befores struc factor calc begins
    rampVal = 1
    rampState = False
  
    # Overwrite properties from file
    if os.path.exists('overwriteProps')== True:
        inp = open('overwriteProps','r')
        ln = inp.readlines(300)
        print(ln)
        for i in range(len(ln)):
            propInd = ln[i].find(' ')
            propName = ln[i][0:propInd]
            propValueStr = ln[i][propInd+3:len(ln[i])-1]  #Property value is last part of line, but drop ' _ ' and the \n
            # check to see if property is one that gets updated here
            #spring_constant, chiN, chiAB, bjerrum_length
            if propName == 'spring_constant' or propName == 'k':
                spring_constant = utl.getValues(propValueStr)
                print('overwriting spring_constant to ', spring_constant)
            elif propName == 'chiN':
                chiN = utl.getValues(propValueStr)
                chiAB = chiN/np.sum(N)
                print('overwriting chiN to ', chiN)
            elif propName == 'chiAB':
                chiAB = utl.getValues(propValueStr)
                print('overwriting chiAB to ', chiAB)
            elif propName == 'bjerrum_length':
                bjerrum_length = utl.getValues(propValueStr)
                print('overwriting bjerrum_length to ', bjerrum_length)
            elif propName == 'delt':
                delt= utl.getValues(propValueStr)
                print('overwriting delta t to ', delt)
            elif propName == 'q_smear':
                smearing_size = utl.getValues(propValueStr)
                print('overwriting smearing size to ', smearing_size)
            elif propName == 'ramp':
                rampVal = utl.getValues(propValueStr)
                rampState = True
                print('Ramping chi by ', rampVal)
            elif propName == 'nodeCount':
                for j in range(Dim):
                    nodeCount[j] = np.floor(utl.getValues(propValueStr))
                print('Overwriting node count to ', nodeCount)
            elif propName == 'nodeRatio':
                nodeRatio = utl.getValues(propValueStr)
                nodeCount = np.floor(nodeRatio*box)
                print('Overwriting node count to ', nodeCount)
                for i1 in range(len(nodeCount)):
                    if (nodeCount[i1] % 2) == 0:
                        nodeCount[i1] += 1
        inp.close()

    # Derived/defined properties from main
    if np.max(saltConc)>0:
        addedSalt = 1
    else:
        addedSalt = 0

    # see if simple fluid
    simpleFluid = False
    if N[0] == [1]: 
        simpleFluid = True
        




    a = box[0]
    otp = open( propertiesName , 'w' )  # open the document to write
    
    #Write Dimensionality
    ln = 'Dim %d\n\n' % Dim
    otp.write( ln )
    
    # Write max_steps
    ln = 'max_steps %d\n' % max_steps
    otp.write( ln )
    
    # Write log frequency
    ln = 'log_freq %d\n' % log_freq
    otp.write( ln )

    # Write binary_freq
    ln = 'binary_freq %d\n' % binary_freq
    otp.write( ln )

    # Write trajectory frequency
    ln = 'traj_freq %d\n' % traj_freq
    otp.write( ln )

    # Write charge info
    if do_charge == 1:
        ln = 'charges %f %f\n\n' % (bjerrum_length, smearing_size)
        otp.write( ln )

    # Write pme order
    ln = 'pmeorder %d\n\n' % pmeorder
    otp.write( ln )

    # Write delta t
    ln = 'delt %f\n' % delt
    otp.write( ln )


    # Write file location for molecule data
    ln = 'read_data input.data\n' #Note- this does not automatically update if var name changes- needed for autogen script
    otp.write( ln )
    # Write structure factor freqency
    if struc_freq != 0:
        if skip_steps != 0:
            ln = 'compute avg_sk 1 freq %d wait %d\n\n' % (struc_freq, skip_steps)
        else:
            ln = 'compute avg_sk 1 freq %d\n\n' % (struc_freq)
        otp.write( ln )
        # Write structure factor freqency
    else:
        ln = '\n'
        otp.write( ln )


    # Write Integrator
    otp.write( 'integrator ' + integratorGroupName + ' ' + integrator + '\n\n')

    # Write resume location if not present
    if readResumeLocation != 'none':
        ln = 'read_resume ' + readResumeLocation + '\n\n'
        otp.write( ln )

    # Write box dimensions
    if Dim > 0:
        ln = 'Nx %d\n' % nodeCount[0]
        otp.write( ln )
    if Dim > 1:
        ln = 'Ny %d\n' % nodeCount[1]
        otp.write( ln )
    if Dim > 2:
        ln = 'Nz %d\n' % nodeCount[2]
        otp.write( ln )

    otp.write('\n')


    # _______Write Bond information__________________
    # all polymers have bond type 1- this is the backbone bond
    if simpleFluid:
        ln = 'bond 1 harmonic %f 0.0\n\n' % spring_constant
        otp.write( ln )
    else:
        temp = 1.5
        ln = 'bond 1 harmonic %f 0.0\n' % temp
        otp.write( ln )
        if np.sum(drudeOscillator)>=1:       #Check if there is a drude-oscillator
            # Add in drude oscillator
            ln = 'bond 2 harmonic %f 0.0\n\n' % spring_constant
            otp.write( ln )
            #if isIonomer == 1:          #Check if there is an ionomer. Since all of our ionomers need a drude-oscillator, we can nest here
                #ln = 'bond 3 harmonic 0.0 0.0\n'
                #otp.write( ln )
        otp.write ('\n')
        # Bonded interactions

    # Establish the number of bonds
    numTypes = np.max (tp)
    # numNonbonded = int(numTypes * (numTypes + 1)/2-numTypes* isIonomer)
    numNonbonded = int(numTypes * (numTypes + 1)/2-numTypes* np.max(drudeOscillator))
    if simpleFluid: numNonbonded = 3; #Too lazy to change my math to account for this
    if addedSalt == 1:
        # add salt nonbonded if there is any
        numNonbonded += 2*numTypes + 3
    ln ='n_gaussians %d\n' % numNonbonded
    otp.write( ln )


    for q1 in range (numTypes):
        for q2 in range (q1, numTypes):             # Loop through all bond pairs
            if (q1 == q2 and q1 + 1 != drude_type):  #Self pair, not drude
                temp = kappa / (2 * rho0)
                ln = 'gaussian %d %d %f %f\n' % (q1 + 1, q2+1, temp, std_dev)
                otp.write( ln )
            elif (q1 != q2 and (q1+1 == numTypes or q2+ 1 == numTypes) and (q1 +1 != drude_type and q2+1 !=drude_type)): #Counterion-AB
                temp =(kappa + chiCounterion)/ rho0
                ln = 'gaussian %d %d %f %f\n' % (q1 + 1, q2+1, temp, std_dev)
                otp.write( ln )

            elif (q1 != q2 and q1 +1 != drude_type and q2 + 1 != drude_type): #A-B interaction
                temp = (kappa + chiAB)/ rho0
                ln = 'gaussian %d %d %f %f' % (q1 + 1, q2+1, temp, std_dev)
                if rampState:
                    temp = (kappa + chiAB*rampVal)/ rho0
                    ln = ln + ' ramp %f\n' % temp
                else:
                    ln = ln +'\n'
                otp.write( ln )
    # Add salt non-bonded
    if addedSalt == 1:
        temp = kappa / (2 * rho0)
        for q1 in range(numTypes):
            # Interactions with all other atoms
            ln = 'gaussian %d %d %f %f\n' % (q1 + 1, numTypes+1, temp, std_dev)
            otp.write( ln )
            ln = 'gaussian %d %d %f %f\n' % (q1 + 1, numTypes+2, temp, std_dev)
            otp.write( ln )
        # Self interactions
        ln = 'gaussian %d %d %f %f\n' % (numTypes+1, numTypes+1, temp, std_dev)
        otp.write( ln )
        ln = 'gaussian %d %d %f %f\n' % (numTypes+2, numTypes+2, temp, std_dev)
        otp.write( ln )
        # + - interactions
        ln = 'gaussian %d %d %f %f\n' % (numTypes+1, numTypes+2, temp, std_dev)
        otp.write( ln )
   
    otp.close()



