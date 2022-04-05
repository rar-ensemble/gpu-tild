#!/home/rrig/Install/anaconda3/bin/python3
import numpy as np
import os, sys
from numpy.random import rand, normal
import polymers as poly
import write_data as wd


# Uses 'b' for dimensionless units

Dim = 2                     # dimensionality of the system
poly_types = 2              # Number of different polymer types
rho0 = 3.0                  # density of the polymer

# Number of monomers of each type in block segments
Np = 15
Nnp = 15

# Properties of each polymer type. Length of outer list must match poly_types and ??? Lengths of inner lists must match each other
tp = 0.95238
phi = np.array([.35,0.65])
q = [[0,0],[0]]
N = [[5,5],[10]]
types = [[1,2],[3]]
wlc_flag = [[1,0],[0]]
drude_q = [[0,0],[0]]
#q = [[0.,0.],[1.],[-1.]]
#N = [[Np, Nnp],[1],[1]]
#types = [[1,2],[3],[3]]                 # Also must be list of lists #1, 2,3 is a b or c (generic- as many as needed)
                                # Sizes of lists must match N, q as it assigns monomer type for each block.
                                # Probably best/should match N exactly (1=Na, 2=Nb...)
#ionomerCondition = [[0,0],[0],[0]]  # If material is an ionomer, and if so, set 1 for which is the charged material
#wlc_flag = [[0,0],[0],[0]]              # Worm-like chain flag
#drude_q = [[0,0],[0],[0]]               # mag of drude charges. 0 disables drude osc on this block

box = np.array([256, 256, 5])                 # Size of the simulation box

Rg = ((Np-1.)/6.)**0.5      # Radius of Gyration
L = 15.0*Rg                 # ???
name = "input.data"         # name of the file with input parameters



#### NOTES ####
# Example i ran- N = [[Na,Nb], [Nashort], [Nbshort]] 2d-bue
# fitst bracket- list of evything
# second bbracke- list of each block lenght- so [[Na,Nb]] gives diblock because second [] has Na and Nb
# [[Na,Nb]] gives a diblock. 

#N =     [[Np, Nnp, Np, Nnp, Np, Nnp, Np, Nnp, Np, Nnp]]   # Must be a list of lists





########################
# END OF INPUT SECTION #
########################

##### Check to see if inputs were configured correctly
if (len(phi) != poly_types):
    sys.exit("Phi configured incorrectly")

if (len(q) != poly_types):
    sys.exit("Charges (q) configured incorrectly")

if (len(drude_q) != poly_types):
    sys.exit("Charges (q) configured incorrectly")

if (len(N) != poly_types):
    sys.exit("Block Lengths (N) configured incorrectly")

# Check If there is a charge- set flag initially to false
do_charge = 0


if (Dim == 2):
    box[2] = 5.0;

# Calculate the volume of the simulation box
V = 1.0
for j in range(0,Dim):
    V = V * box[j]

# ??? Storage for number of molecules and blocks
nmolecs = np.zeros(poly_types,'i')
nblocks = np.zeros(poly_types,'i')
ndrude_sites = np.zeros(poly_types,'i')


nstot = 0                       # Total number of monomer units
nbonds_tot = 0                  # Total number of monomer bonds
nangles_tot = 0
max_type = -3

for i in range(poly_types):
    # First, figure out if we're dealing with a block polymer. This does this by looking at the length of each secondary
    # list in N. This is set up so it gives the length of each blocks in that Monomer Type
    nblocks[i] = len(N[i])

    # Check for a charge. If there is, set flag do_charge to true
    nangs_molec = 0;

    ndrude_sites[i] = 0;
    for j in range(nblocks[i]):
        # only enter loop if charge condition is false
        if ( do_charge == 0.):             
            if ( q[i][j] != 0. ):
                do_charge = 1
            
        # Accumulate drude_sites as needed
        if ( drude_q[i][j] != 0. ):
            ndrude_sites[i] = ndrude_sites[i] + N[i][j]
            if ( do_charge == 0 ):
              do_charge = 1;

                
        if ( wlc_flag[i][j] == 1 ):
            if ( j < nblocks[i]-1 and wlc_flag[i][j+1] == 1 ):
                nangs_molec = nangs_molec + N[i][j]
            else:
                nangs_molec = nangs_molec + N[i][j] - 2

        if ( types[i][j] > max_type ):
            max_type = types[i][j];
    Ntot = np.sum(N[i])             # Overal length of poly_type i
    
    # Generate total number of molecules based on volume, density, and length of each
    nmolecs[i] = int( phi[i] * rho0 * V / Ntot )


    # updaterunning count of total number of monomer units and bonds
    nstot = nstot + int( nmolecs[i] * (Ntot + ndrude_sites[i]) )
    nbonds_tot = nbonds_tot + int( nmolecs[i] * (Ntot-1+ndrude_sites[i]) )
    nangles_tot = nangles_tot + nangs_molec * nmolecs[i]

    # Add drude bonds, sites to the running tally


print('# Molecules of each Type:', nmolecs)
print('Total particles:', nstot, 'Total bonds:', nbonds_tot, 'Total angles:', nangles_tot)
if ( do_charge ):
    print("Implementing charges!!!");

# Initialize storage variables        
x = np.zeros((nstot,3),'d')
bonds = np.zeros((nbonds_tot,3),'i')
angles = np.zeros((nangles_tot,4),'i')
mID = np.zeros(nstot,'i')
typ = np.zeros(nstot,'i')
ind = np.zeros(nstot,'i')
charges = np.zeros(nstot,'d')

bind = 0
aind = 0
ind_shift = 0
mol_ind = 0
for i in range(poly_types):
    
    Ntot = np.sum(N[i])         # Total length of the polymer
    Nlist = N[i]                # Length of each block
    nblocks = len(Nlist)        # Number of blocks
    qlist = q[i];               # Charge of each blocka
    loc_drude_q = drude_q[i]
    loc_types = types[i]
    wlc_list = wlc_flag[i]

    for j in range(nblocks):
        if ( loc_drude_q[j] != 0. ):
            loc_types.append( max_type+1 )
            break;


    print("First molecule type has", nblocks, " blocks.")
        
    for j in range(nmolecs[i]):
        #################################################################################
        # IONOMER PLAN- Introduce a variable that calls the ionomer building function when I set flags (to be added) before
        # Another possible plan/ way to introduce- use the same generation, then add the counterion after with flag
        # call make_diblock from polymers.py

        xl, bl, angl, tp, ql = poly.make_nblock(Nlist, loc_types, qlist, loc_drude_q, wlc_list, box, Dim)

        mol = np.zeros(Ntot+ndrude_sites[i],'i')                #what ? polymer number it is
        mol[:] = mol_ind + 1
        mol_ind += 1                            # increment polymer number
        
        il = np.zeros(Ntot+ndrude_sites[i],'i')                 # Global
        for k in range(Ntot+ndrude_sites[i]):
            il[k] = ind_shift + k + 1           # Monomer number in polymer- global with ind_shift from local number from make_diblock
        

        bl[:,1] = bl[:,1] + ind_shift           # Shift bond info to global index number
        bl[:,2] = bl[:,2] + ind_shift           # Shift bond info to global index number

        for k in range(Ntot-1+ndrude_sites[i]):
            bonds[bind,:] = bl[k,:]             # global bond storage
            bind = bind + 1                     # global bond storage index

        n_new_angs = np.shape(angl)[0]
        if ( n_new_angs > 0 ):
            angl[:,1] = angl[:,1] + ind_shift
            angl[:,2] = angl[:,2] + ind_shift
            angl[:,3] = angl[:,3] + ind_shift
            #print(angl)

            for k in range(n_new_angs):
                angles[aind,:] = angl[k,:]
                aind = aind + 1

        x[ind_shift:ind_shift+Ntot+ndrude_sites[i],:] = xl      # global storage of position
        mID[ind_shift:ind_shift+Ntot+ndrude_sites[i]] = mol     # global storage of ? polymer number
        typ[ind_shift:ind_shift+Ntot+ndrude_sites[i]] = tp      # global storage of type number
        ind[ind_shift:ind_shift+Ntot+ndrude_sites[i]] = il      # global storage of monomer number
        charges[ind_shift:ind_shift+Ntot+ndrude_sites[i]] = ql  # global storage of charge
                    
        ind_shift += Ntot+ndrude_sites[i]                # index shift for global count of monomers


if ( np.sum(charges) != 0.0 ):                  # Check to make sure
    print("WARNING WARNING!!!\n\nWARNING!! WARNING!!\nSystem not charge neutral!! Q = ", np.sum(charges));

wd.write_data(name, ind, mID, typ, charges, box, x, bonds, angles)  #Write the data to the file

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
