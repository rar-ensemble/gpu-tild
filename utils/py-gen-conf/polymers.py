#!/home/rrig/Install/anaconda3/bin/python3
import numpy as np
import os, sys
from numpy.random import rand, normal


def make_nblock(Nlist, types, qlist, drude_q, wlc_flag, box, Dim):
    # What I think this does is generate a block polymer (not necissarily diblock)s
    Ntot = np.sum(Nlist)            # Total length of molecule (including all blocks)
    nblocks = len(Nlist)            # Number of blocks

    if ( len(types) != len(Nlist) ):
      print("HORRIBLE PROBLEM!!!")
      exit(1)
    
    ndrude = 0;
    for i in range(nblocks):
        if ( drude_q[i] != 0. ):
            ndrude += Nlist[i]


    nbond_types = 1;
    for i in range(nblocks-1):
        if ( wlc_flag[i] != wlc_flag[i+1]):
            nbond_types = 2;

    # Storage arrays
    typ = np.zeros(Ntot+ndrude,'i')            # Monomer type for each monomer
    r = np.zeros((Ntot+ndrude, 3),'d')          # Position vectory for each monomer
    bonds = np.zeros((Ntot-1+ndrude, 3), 'i')    # What monomers are bonded (relative number for monomer) (one less bond that monomers)
    charges = np.zeros(Ntot+ndrude, 'd')        # Charge of the monomer

    nangles = 0
    for i in range(nblocks):
        if ( wlc_flag[i] == 1 ):
            if ( i == nblocks-1 or wlc_flag[i+1] == 0 ):
                nangles += Nlist[i]-2;
            else:
                nangles += Nlist[i]
#    print("Molecule has", nangles, "angles")
    angles = np.zeros((nangles,4),'i')

    ind = 0;
    for j in range(0,Dim):
        r[ind,j] = rand() * box[j]        # Random starting position for each monomer

    typ[0] = types[0]                   # First monomer is the first type (guaranteed)
    charges[0] = qlist[0]               # First monomer has charge of first type

    ind = 1;


    bond_ind = 0;
    if ( Ntot > 1 ):                    # If there is more than one monomer in polymer, assign bond info for first bond
      if ( nbond_types == 2 and wlc_flag[0] == 1 ):
        bonds[0,0] = 2
      else:
        bonds[0,0] = 1
      bonds[0,1] = 1;                   # monomer 1 in bond
      bonds[0,2] = 2;                   # monomer 2 in bond
      bond_ind = 1

    block_ind = 1                       # Index of what block is being built. Used to assign type info

    ang_ind = 0
    if ( wlc_flag[0] == 1 and Nlist[0] >= 3 ):
        angles[0,0] = 1
        angles[0,1] = 1
        angles[0,2] = 2
        angles[0,3] = 3
        ang_ind = 1



    for i in range(1,Ntot):
        for j in range(0,Dim):
            r[ind,j] = r[ind-1,j] + normal();   # next mono position is some distance away from previous mono position
            if ( r[ind,j] > box[j] ):
                r[ind,j] = r[ind,j] - box[j]    # check PBC
            elif (r[ind,j] < 0.0 ):
                r[ind,j] = r[ind,j] + box[j]    # check PBC

        if ( i >= np.sum(Nlist[0:block_ind] ) ):
            block_ind += 1                  # if i is passed current block length, advance block index 1

        typ[ind] = types[block_ind-1]         # assign monomer to the type for the block
        charges[ind] = qlist[block_ind-1]        # assign charge for the monomer type
        ind += 1

        if ( i < (Ntot-1)):                 # Generate bond info for bond starting at this bond
            #bonds[bond_ind,0] = 1;                 # bond type
            if ( nbond_types == 2 and wlc_flag[block_ind-1] == 1 ):
              bonds[bond_ind,0] = 2
            else: 
              bonds[bond_ind,0] = 1

            bonds[bond_ind,1] = i+1;               # monomer 1 in bond
            bonds[bond_ind,2] = i+2;               # monomer 2 in bond
            bond_ind += 1

        if ( wlc_flag[block_ind-1] == 1  and \
                ( i < np.sum(Nlist[0:block_ind])-2 or \
                ( block_ind < nblocks and wlc_flag[block_ind] == 1 ) ) ):
                
            #print(i, np.sum(Nlist[0:block_ind])-2)
            #print( i< np.sum(Nlist[0:(block_ind)])-2,  block_ind < nblocks, wlc_flag[block_ind-1] == 1)
            angles[ang_ind,0] = 1
            angles[ang_ind,1] = i+1
            angles[ang_ind,2] = i+2
            angles[ang_ind,3] = i+3
            ang_ind = ang_ind + 1

    # Add the drude particles
    Nstart = 0;
    for i in range(nblocks):
        if ( drude_q[i] != 0. ):
            for k in range(Nstart, Nstart+Nlist[i]):
                for j in range(Dim):
                    r[ind,j] = r[k,j] + normal();
                    if ( r[ind,j] > box[j] ):
                        r[ind,j] = r[ind,j] - box[j]    
                    elif (r[ind,j] < 0.0 ):
                        r[ind,j] = r[ind,j] + box[j]  

                typ[ind] = types[nblocks];
                charges[k] = drude_q[i]
                charges[ind] = -drude_q[i];
                ind += 1

                bonds[bond_ind,0] = 2;
                bonds[bond_ind,1] = k;
                bonds[bond_ind,2] = ind;
                bond_ind += 1  

    

    return r, bonds, angles, typ, charges
