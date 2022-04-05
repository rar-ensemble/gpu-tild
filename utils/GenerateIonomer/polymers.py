#!/home/rrig/Install/anaconda3/bin/python3
import numpy as np
import os, sys
from numpy.random import rand, normal
import utilities


def make_nblock(Nlist, types, qlist, box, Dim, ionomerCondition, ionomerCharge,drudeOscillator):
    # What I think this does is generate a block polymer (not necissarily diblock)s
    Ntot = np.sum(Nlist)            # Total length of molecule (including all blocks)
    nblocks = len(Nlist)            # Number of blocks
    
    if ( len(types) != len(Nlist) ):
      print("HORRIBLE PROBLEM!!!")
      exit(1)
    # ionomer counterion

    # Storage arrays
    typ = np.zeros(Ntot,'i')            # Monomer type for each monomer
    r = np.zeros((Ntot,3),'d')          # Position vectory for each monomer
    bonds = np.zeros((Ntot-1,3),'i')    # What monomers are bonded (relative number for monomer) (one less bond that monomers)
    charges = np.zeros(Ntot,'d')        # ??? Charge of the monomer

    for j in range(0,Dim):
        r[0,j] = rand() * box[j]        # Random starting position for each monomer

    typ[0] = types[0]                   # First monomer is the first type (guaranteed)

    if ( Ntot > 1 ):                            # If there is more than one monomer in polymer, assign bond info for first bond
      bonds[0,0] = 1 + drudeOscillator[0];      # bond type
      bonds[0,1] = 1;                           # monomer 1 in bond
      bonds[0,2] = 2;                           # monomer 2 in bond
      charges[0] = qlist[0]*((-1)**(ionomerCondition[0]+1))   # Formula to cancluate charges
      simpleFluid = False   # Flag for simple fluid
    if (Ntot ==2):  #condition for simple fluids
      bonds[0,1] = 1;
      for j in range(0,Dim):
          r[1,j] = utilities.periodicBC1(r[0,j] + normal(),box[j]) # next mono position is some distance away from mono 2, checked for PBC
      charges[1] = -1*charges[0]
      typ[1] = 2
      simpleFluid = True

    block_ind = 1                       # Index of what block is being build. Used to assign type info
    for i in range(1,Ntot):    #I have no idea why it is range(1,Ntot), but it seems to work fine...
        if simpleFluid: break   # Exit out if simple fluid
        if (drudeOscillator[block_ind-1]==1): # This section runs if there is a drude-oscillator
            # This checks to see if it is a charged ionomer section. If so, every other monomer placed is dipole and last is counterion
            placementCounter = np.sum(Nlist[0:block_ind])-i     # keep track of what is being placed in iono
            if (i==np.sum(Nlist[0:block_ind])-1 and ionomerCondition[block_ind-1]==1):               # placing counterion for backbone
                for j in range(0,Dim):
                    r[i,j] = r[i-np.sum(Nlist[0:block_ind])+3,j] + normal();   # next mono position is some distance away from mono 2
                    r[i,j] = utilities.periodicBC1(r[i,j],box[j])   # Check PBC
                    
                typ[i] = np.max(types) + 2          # counter ion has type two above max input type number
                charges[i] = 1 * ionomerCharge        # equal in mag to all of the other charges, opp sign
                charges[i-Nlist[block_ind-1]+3] -= ionomerCharge  # Add counterion charge
                bonds[i,0] = 3                      # different bond type than backbone
                bonds[i,1] = i - Nlist[block_ind-1]+ 4
                bonds[i,2] = i + 1
                block_ind += 1          # advance the block index since we know this is the last one
            else:
                # Since there is always an odd number for placement counter (2x monomers for dipole and 1x counterion), we can set
                # odd counters to backbone and even for the dipole in Drude-oscillator model. I am checking even/odd by raising -1 to the power in placment
                # if it is an even counter, then it places it next to the previous as the dipole
                # if it is an odd counter, then it places it next two 2 back (previous backbone monomer)
                placeCheck =int( 0.5*(1-(-1)**(placementCounter+1+ionomerCondition[block_ind-1])))      # Odd = 1, even = 0
                for j in range(0,Dim):
                    r[i,j] = r[i-1-placeCheck,j] + normal();   # next mono position is some distance away from mono 2
                    r[i,j] = utilities.periodicBC1(r[i,j],box[j])   # Check PBC

                typ[i] = np.max(types)*(1-placeCheck) + 1          # drude-oscillator
                charges[i] = ((-1)**(placementCounter+1))*qlist[block_ind-1]         # equal in mag to all of the other charges, opp sign
                #print(charges)
                bonds[i,0] = 1 + placeCheck                      # different bond type than backbone
                bonds[i,1] = i + placeCheck
                if (placementCounter != 2 and ionomerCondition[block_ind-1]==1):
                    bonds[i,2] = i + 2
                elif ionomerCondition[block_ind-1]==1:
                    bonds[i,2] = i + 3                 # Account for counterion
                elif placementCounter != 1 and ionomerCondition[block_ind-1]==0:      #advance block index for non-ionomers
                    bonds[i,2] = i + 2
                else:                       # Corrections for non ionomer DO
                    bonds[i,2] = i + 2
                    bonds[i,0] = 1
                    #typ[i] = types[block_ind]
                    #charges[i] =qlist[block_ind]
                    #print(charges)
                    block_ind+=1
                    
        else:                                         # non-drude oscillator generation
            for j in range(0,Dim):
                r[i,j] = r[i-1,j] + normal();   # next mono position is some distance away from previous mono position
                r[i,j] = utilities.periodicBC1(r[i,j],box[j])   # Check PBC
            a = np.sum(Nlist[0:block_ind] )
            
            #if ( i >= np.sum(Nlist[0:block_ind] ) ):
                #block_ind += 1                  # if i is passed current block length, advance block index 1

            typ[i] = types[block_ind-1]         # assign monomer to the type for the block
            charges[i] = qlist[block_ind-1]        # assign charge for the monomer type #### Added [i]
            #print(charges)
            if ( i < (Ntot-1)):                 # Generate bond info for bond starting at this bond
                if (ionomerCondition[block_ind-1]==1 and i==np.sum(Nlist[0:block_ind])-2):
                    bonds[i,0] = 1;                 # bond type
                    bonds[i,1] = i+1;               # monomer 1 in bond
                    bonds[i,2] = i+3;
                else:
                    bonds[i,0] = 1 +ionomerCondition[block_ind-1];                 # bond type
                    bonds[i,1] = i+1;               # monomer 2 in bond
                    bonds[i,2] = i+2;               # monomer 2 in bond
            a = np.sum(Nlist[0:block_ind] )

            if ( i >= np.sum(Nlist[0:block_ind] )-1 ):
                block_ind += 1                  # if i is passed current block length, advance block index 1
    #Remove charge bond
    bonds = np.delete(bonds,np.where(bonds[:,0]==3),axis=0)

    if np.unique(bonds[:,0].shape[0])==1:   #Not tested complex
        bonds[:,0]=1

    #print(charges)
    #print(bonds)
    #print(typ)
    return r, bonds, typ, charges

