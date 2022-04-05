#!/home/rrig/Install/anaconda3/bin/python3
import numpy as np
import os, sys
from numpy.random import rand, normal
import polymers as poly
import write_data as wd
import write_prop as wp
import utilities as utl

## NOTE- this is modified main script to be run from generate_inputs
# The difference is this is no longer what is called, all of it is in def, and there is an overwrite script
# Overwrite script seaches for variables in file created by gen_input that need to be overwritten for mass generation
# of different system properties

# Uses 'b' for dimensionless units
def genSystem():
    Dim = 3                     # dimensionality of the system
    poly_types = 1              # Number of different polymer types
    rho0 = 3.85                 # density of the polymer

    # Number of monomers of each type in block segments
    Na = 5
    Nb = 15

    # Properties of each polymer type. Length of outer list must match poly_types and ??? Lengths of inner lists must match each other
    phi = np.array([1])           # Phi for each poly_type (volume fraction)
    q = [[-1,0,-1,0]]               # Charge of each polymer by block. Used for dipole charge in ionomer of backbone
    N = [[Na,Nb,Na,Nb]]          # Must be a list of lists. Secondary lists are lengths of each block for that polymer type
    types =[[1,2,1,2]]            # Also must be list of lists #1, 2,3 is a b or c (generic- as many as needed)
                                    # Sizes of lists must match N, q as it assigns monomer type for each block.
                                    # Probably best/should match N exactly (1=Na, 2=Nb...)
                                    # Code assumes sequential (1,2,...) as inomer counterion is next in list

    drudeOscillator = [[1,0,1,0]]   # Flag to add drude-oscillator to a polymer. Uses q for charges to net 0

    ionomerCondition = [[1,0,1,0]] # If material is an ionomer, and if so, set 1 for which is the charged material
                                    # Note that it needs to be either 1 or 0 for code to work. Set it to be all zeros for no ionomer
    ionomerCharge = 1              # Charge for ionomer/counterion

    relativeSaltConc = [0] # Add salt proportional to number of molecules of a type: Na, Nb, Nc... or types [1,2,3,...]
                                # code works off first number being proportional to type 1, then type 2
    saltCharge = -1              # Charge of the salt

    box = np.array([40, 40,40])                 # Size of the simulation box

    Rg = ((Na-1.)/6.)**0.5      # Radius of Gyration
    L = 15.0*Rg                 # ???
    name = 'input.data'         # name of the file with molecule data
    propertiesName = 'input'    # name for system properties

    #### NOTES ####
    # Example i ran- N = [[Na,Nb], [Nashort], [Nbshort]] 2d-bue
    # fitst bracket- list of evything
    # second bbracke- list of each block lenght- so [[Na,Nb]] gives diblock because second [] has Na and Nb
    # [[Na,Nb,Na,Nb,Na,Nb]] gives a multiblock of repeating a and b 

    #N =     [[Np, Nnp, Np, Nnp, Np, Nnp, Np, Nnp, Np, Nnp]]   # Must be a list of lists
    
    #check for any overwrites
    if os.path.exists('overwriteProps')== True:
        inp = open('overwriteProps','r')
        ln = inp.readlines(100)
        print(ln)
        for i in range(len(ln)):
            propInd = ln[i].find(' ')
            propName = ln[i][0:propInd]
            propValueStr = ln[i][propInd+3:len(ln[i])-1]  #Property value is last part of line, but drop ' _ ' and the \n
            # check to see if property is one that gets updated here
            if propName == 'direc':  #This one always should be there- directory to create the files
                print('writing in directory', propValueStr)
                name = propValueStr +'/' + name
                propertiesName = propValueStr + '/' + propertiesName
                print(name)
                print(propertiesName)
            elif propName == 'rho0':
                rho0 = utl.getValues(propValueStr)
                print('overwriting rho0 to ', rho0)
            elif propName == 'Na':
                Na = int(utl.getValues(propValueStr))
                N = [[Na,Nb,Na,Nb]]
                print('overwriting Na to ', Na)
            elif propName == 'Nb':
                Nb = int(utl.getValues(propValueStr))
                N = [[Na,Nb,Na,Nb]]
                print('overwriting Nb to ', Nb)
            elif propName == 'box' or propName == 'vol':
                box[box>0] = utl.getValues(propValueStr)
                print('overwriting box to ', box)
            elif propName == 'saltA':
                relativeSaltConc[0] = utl.getValues(propValueStr)
                print('overwriting relative salt concentraion of atom type 1 to ', relativeSaltConc)
            elif propName == 'saltB':
                relativeSaltConc[1] = utl.getValues(propValueStr)
                print('overwriting relative salt concentraion of atom type 2 to ', relativeSaltConc)
            elif propName == 'qA':
                q[0][0] = utl.getValues(propValueStr)   #ThIS DOES NOT WORK GENERALLY
                print('overwriting charge of A to ', utl.getValues(propValueStr))
            elif propName == 'frac1A':
                Nb = int(Na/ utl.getValues(propValueStr))
                N = [[Na,Nb,Na,Nb]]
                print('overwriting number of B to ', Nb)

        inp.close()

    ########################
    # END OF INPUT SECTION #
    ########################

    ##### Check to see if inputs were configured correctly
    if (len(phi) != poly_types):
        sys.exit("Phi configured incorrectly")

    if (len(q) != poly_types):
        sys.exit("Charges (q) configured incorrectly")

    if (len(N) != poly_types):
        sys.exit("Block Lengths (N) configured incorrectly")

    if (np.sum(phi)!= 1):
        sys.exit("Volume fractions (phi) do not add to zero")

    # Check If there is a charge- set flag initially to false
    do_charge = 0

    # Check if there is an ionomer
    isIonomer = 0
    if (np.max(ionomerCondition)>0):
        isIonomer=1             # flag to build ionomer

    if (Dim == 2):
        box[2] = 5.0;

    # Calculate the volume of the simulation box
    V = 1.0
    for j in range(0,Dim):
        V = V * box[j]

    # ??? Storage for number of molecules and blocks
    nmolecs = np.zeros(poly_types,'i')
    nblocks = np.zeros(poly_types,'i')


    nstot = 0                       # Total number of monomer units
    nbonds_tot = 0                  # Total number of monomer bonds

    for i in range(poly_types):
        # First, figure out if we're dealing with a block polymer. This does this by looking at the length of each secondary
        # list in N. This is set up so it gives the length of each blocks in that Monomer Type
        nblocks[i] = len(N[i])

        # Check for a charge. If there is, set flag do_charge to true
        for j in range(nblocks[i]):
            # only enter loop if charge condition is false
            if ( do_charge == 0.):             
                if ( q[i][j] != 0. ):
                    do_charge = 1
                    # exit loop once charge is determined to be presesent
                    break 
        #ionomerAddition=np.add(ionomerCondition[i], np.prod([N[i],ionomerCondition[i]],axis=0))
        additionalComponents = np.add(ionomerCondition[i],np.prod([N[i], drudeOscillator[i]],axis = 0))    # Additional components beyond base bonomer
        Ntot = np.sum(N[i])  # Overal length of poly_type i including ionomers without additions       
    
        # Generate total number of molecules based on volume, density, and length of each
        nmolecs[i] = int( phi[i] * rho0 * V / Ntot )

        # Update Ntot to include additional components
        Ntot += np.sum(additionalComponents)

        # updaterunning count of total number of monomer units and bonds
        nstot = nstot + int( nmolecs[i] * (Ntot))
        nbonds_tot = nbonds_tot + int( nmolecs[i] * (Ntot-1-np.sum(ionomerCondition)) )                     # Note that Ntot acounts for ionomer inclusion
                                                                                    # -1 as there is one less bond than molecules
                                                                                    # - sum IonCond drops counterion bond
        # This includes for the counterion in the data storage if there is an ionomer. 
        # Multiplies the number of ionomer sections for the polymer tupe (sum)
        # By the number of monomers of that type (nmolecs[i])

    # if there is an ionomer, update nstot to include counterion

    print('# Molecules of each Type:', nmolecs)
    print('Total particles:', nstot, 'Total bonds:', nbonds_tot)
    if ( do_charge ):
        print("Implementing charges!!!");

    # Initialize storage variables        
    x = np.zeros((nstot,3),'d')
    bonds = np.zeros((nbonds_tot,3),'i')
    mID = np.zeros(nstot,'i')
    typ = np.zeros(nstot,'i')
    ind = np.zeros(nstot,'i')
    charges = np.zeros(nstot,'d')

    bind = 0
    ind_shift = 0
    mol_ind = 0

    for i in range(poly_types):
        additionalComponents = np.add(ionomerCondition[i],np.prod([N[i], drudeOscillator[i]],axis = 0))
        Ntot = np.sum(N[i]) + np.sum(additionalComponents)        # Total length of the polymer
        Nlist = np.add(N[i], additionalComponents)             # Length of each block w/ ionomer counterion added
        nblocks = len(Nlist)        # Number of blocks
        qlist = q[i];               # Charge of each block

        print("First molecule type has", nblocks, " blocks.") 
        for j in range(nmolecs[i]):
            #################################################################################
            # IONOMER PLAN- Introduce a variable that calls the ionomer building function when I set flags (to be added) before
            # Another possible plan/ way to introduce- use the same generation, then add the counterion after with flag
            # call make_diblock from polymers.py
            xl, bl, tp, ql = poly.make_nblock(Nlist, types[i], qlist, box, Dim, ionomerCondition[i],ionomerCharge,drudeOscillator[i])
            mol = np.zeros(Ntot,'i')                #what ? polymer number it is
            mol[:] = mol_ind + 1
            mol_ind += 1                            # increment polymer number

            il = np.zeros(Ntot,'i')                 # Global
            for k in range(Ntot):
                il[k] = ind_shift + k + 1           # Monomer number in polymer- global with ind_shift from local number from make_diblock
        

            bl[:,1] = bl[:,1] + ind_shift           # Shift bond info to global index number
            bl[:,2] = bl[:,2] + ind_shift           # Shift bond info to global index number
            if j == 0:
                blSaved = bl
                tpSaved = tp

            for k in range(bl.shape[0]): #range(Ntot-1-np.sum(ionomerCondition)):
                bonds[bind,:] = bl[k,:]             # global bond storage
                bind = bind + 1                     # global bond storage index

            x[ind_shift:ind_shift+Ntot,:] = xl      # global storage of position
            mID[ind_shift:ind_shift+Ntot] = mol     # global storage of ? polymer number
            typ[ind_shift:ind_shift+Ntot] = tp      # global storage of type number
            ind[ind_shift:ind_shift+Ntot] = il      # global storage of monomer number
            charges[ind_shift:ind_shift+Ntot] = ql  # global storage of charge
                    
            ind_shift += Ntot                       # index shift for global count of monomers

    # Section to add the salt
    if np.max(relativeSaltConc)>0:   #only run this code if there is salt to be added
        print('Adding Salt to the system')
        moleculeTypeCount = np.unique(typ)   #Get the number of molecule types in the system
        saltNumber = 0
        for i1 in range(len(relativeSaltConc)):   #For each relative salt conc value
            typeCount = np.count_nonzero(typ == moleculeTypeCount[i1])
            saltNumber += int(typeCount *relativeSaltConc[i1])
        saltIndex = np.max(moleculeTypeCount)+1 #index value
        moleculeIndex = np.max(mID)+ 1
        xlSalt, tpSalt, qlSalt, mIDSalt = utl.addSalt(box,Dim, saltNumber, saltIndex, saltCharge,moleculeIndex)
        #Append salt data to the end of relevant variables
        x = np.append(x,xlSalt,axis=0)             
        typ = np.append(typ,tpSalt)
        charges = np.append(charges, qlSalt)
        mID = np.append(mID,mIDSalt)
        ind = np.append(ind,np.arange(np.max(ind)+1,np.max(ind)+1+saltNumber*2))   #Add onto total count index

    if ( abs(np.sum(charges)) > 0.00000001 ):                  # Check to make sure system is neutral. Abs used for float rounding error
        print("WARNING WARNING!!!\n\nWARNING!! WARNING!!\nSystem not charge neutral!!\n");

    wd.write_data(name, ind, mID, typ, charges, box, x, bonds, [])  #Write the data to the file
    wp.write_prop( propertiesName, Dim, name, rho0, box, blSaved ,tpSaved, isIonomer,drudeOscillator, do_charge,relativeSaltConc,N)

    # Print Some informtation
    print(name,'successfully written. Some parameters that may be of interest, based on N of first polymer:\n')
    Ns = np.sum(N[0])
    Ns = N[0][0];
    print('Using', Ns, 'to derive parameters below')
    Rg = ((Ns-1.)/6.)**0.5

    print('Rg =', Rg)
    print('Box dimensions in Rg units:', box/Rg)

    C = rho0 * Rg**3/Ns
    print('C = rho0*Rg**3/N[0] =', C)


    kappa = 50 / Ns
    chi = np.array([5., 10., 20., 40., 80.]) / Ns
    print('\nkappa*N[0] =',kappa*Ns,', chiAB*N[0] =', chi*Ns, 'use A_ii =',
          kappa/2/rho0, 'and A_ij =', (kappa + chi)/rho0 )

    kappa = 200 / Ns
    print('\nkappa*N[0] =',kappa*Ns,', chiAB*N[0] =', chi*Ns, 'use A_ii =',
          kappa/2/rho0, 'and A_ij =', (kappa + chi)/rho0 )
    return ()

#Run if this script is called
if __name__ == '__main__':
    genSystem()
