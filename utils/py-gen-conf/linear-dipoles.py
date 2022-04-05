#!/home/rrig/Install/anaconda3/bin/python3
import numpy as np
import os, sys
from numpy.random import rand, normal
import polymers as poly
import write_data as wd


# Uses 'b' for dimensionless units

Dim = 3                     # dimensionality of the system
poly_types = 1              # Number of different polymer types
rho0 = 3.0                  # density of the polymer

# Number of monomers of each type in block segments
Na = 5
Nb = 15

# Properties of each polymer type. Length of outer list must match poly_types and ??? Lengths of inner lists must match each other
phi = np.array([0.5])       # Phi for each poly_type
q = [[-1, 0, -1, 0]]               # Charge of each polymer
N = [[Na, Nb, Na, Nb]]      # Must be a list of lists. Secondary lists are lengths of each block for that polymer type
types = [[1, 2, 1, 2]]            # Also must be list of lists #1, 2,3 is a b or c (generic- as many as needed)
                                # Sizes of lists must match N, q as it assigns monomer type for each block.
                                # Probably best/should match N exactly (1=Na, 2=Nb...)
ionomerCondition = [[1,0,1,0]]  # If material is an ionomer, and if so, set 1 for which is the charged material

box = np.array([125, 125, 125])                 # Size of the simulation box

Rg = ((Na-1.)/6.)**0.5      # Radius of Gyration
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
    Ntot = np.sum(N[i])             # Overal length of poly_type i
    
    # Generate total number of molecules based on volume, density, and length of each
    nmolecs[i] = int( phi[i] * rho0 * V / Ntot )

    # updaterunning count of total number of monomer units and bonds
    nstot = nstot + int( nmolecs[i] * Ntot)
    nbonds_tot = nbonds_tot + int( nmolecs[i] * (Ntot-1) )


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
    
    Ntot = np.sum(N[i])         # Total length of the polymer
    Nlist = N[i]                # Length of each block
    nblocks = len(Nlist)        # Number of blocks
    qlist = q[i];               # Charge of each block

    print("First molecule type has", nblocks, " blocks.")
        
    for j in range(nmolecs[i]):
        #################################################################################
        # IONOMER PLAN- Introduce a variable that calls the ionomer building function when I set flags (to be added) before
        # Another possible plan/ way to introduce- use the same generation, then add the counterion after with flag
        # call make_diblock from polymers.py
        xl, bl, tp, ql = poly.make_nblock(Nlist, types[i], qlist, box, Dim)

        mol = np.zeros(Ntot,'i')                #what ? polymer number it is
        mol[:] = mol_ind + 1
        mol_ind += 1                            # increment polymer number
        
        il = np.zeros(Ntot,'i')                 # Global
        for k in range(Ntot):
            il[k] = ind_shift + k + 1           # Monomer number in polymer- global with ind_shift from local number from make_diblock
        

        bl[:,1] = bl[:,1] + ind_shift           # Shift bond info to global index number
        bl[:,2] = bl[:,2] + ind_shift           # Shift bond info to global index number

        for k in range(Ntot-1):
            bonds[bind,:] = bl[k,:]             # global bond storage
            bind = bind + 1                     # global bond storage index

        x[ind_shift:ind_shift+Ntot,:] = xl      # global storage of position
        mID[ind_shift:ind_shift+Ntot] = mol     # global storage of ? polymer number
        typ[ind_shift:ind_shift+Ntot] = tp      # global storage of type number
        ind[ind_shift:ind_shift+Ntot] = il      # global storage of monomer number
        charges[ind_shift:ind_shift+Ntot] = ql  # global storage of charge
                    
        ind_shift += Ntot                       # index shift for global count of monomers


if ( np.sum(charges) != 0.0 ):                  # Check to make sure
    print("WARNING WARNING!!!\n\nWARNING!! WARNING!!\nSystem not charge neutral!!\n");

wd.write_data(name, ind, mID, typ, charges, box, x, bonds, [])  #Write the data to the file

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


"""
nD = 2; #int(rho0 * V * phiD / N)
ns = nD * (Na+Nb)
n_diblock_bonds = nD*(N-1);

Nh = 10
nH = 0; #int((1.0-phiD)*rho0*V/Nh)
nsH = nH*Nh
n_homo_bonds = nH*(Nh-1)

nP = 2
nG = 5
Nga = 5
Ngb = 5
ns_per_gnp = (Nga + Ngb + 1) * nG + 1 
ns_gnp = nP * ns_per_gnp

bonds_per_gnp = nG*(Nga+Ngb)
n_gnp_bonds = nP*bonds_per_gnp
if ( nP == 0 ):
    ns_gnp = 0
Rp = 2.5

print('nD:', nD, nH, nP)
print('Generating config with', ns, 'diblock sites,', nsH,'homopolymer sites,',ns_gnp,'grafted particle sites')

nstot = ns + nsH + ns_gnp;
nbonds_tot = n_diblock_bonds + n_homo_bonds + n_gnp_bonds
print('Total bonds:',nbonds_tot)

x = np.zeros((nstot,3),'d')
bonds = np.zeros((nbonds_tot,3),'i')
mID = np.zeros(nstot,'i')
typ = np.zeros(nstot,'i')
ind = np.zeros(nstot,'i')

bind = 0
ind_shift = 0
for i in range(0,nD):
    xl, bl, tp = poly.make_diblock(Na, Nb, box, Dim)
    
    il = np.zeros(N,'i')
    mol = np.zeros(N,'i')
    
    ind_shift = i * N
    for j in range(N):
        il[j] = ind_shift + j + 1
        mol[j] = i+1

    bl[:,1] = bl[:,1] + ind_shift
    bl[:,2] = bl[:,2] + ind_shift

    for j in range(N-1):
        bonds[bind,:] = bl[j,:]
        bind = bind + 1

    for j in range(0,N):
        x[ind_shift+j,0] = xl[j,0]
        x[ind_shift+j,1] = xl[j,1]
        x[ind_shift+j,2] = xl[j,2]
        mID[ind_shift+j] = mol[j]
        typ[ind_shift+j] = tp[j]
        ind[ind_shift+j] = il[j]


for i in range(0,nH):
    xl, bl, tp = make_diblock(Nh, 0, box, Dim)
    
    il = np.zeros(N,'i')
    mol = np.zeros(N,'i')
    
    ind_shift = nD*N + i*Nh
    

    bl[:,1] = bl[:,1] + ind_shift
    bl[:,2] = bl[:,2] + ind_shift

    for j in range(Nh-1):
        bonds[bind,:] = bl[j,:]
        bind = bind + 1

    for j in range(Nh):
        il[j] = ind_shift + j + 1
        mol[j] = i+1+nD
        x[ind_shift+j,0] = xl[j,0]
        x[ind_shift+j,1] = xl[j,1]
        x[ind_shift+j,2] = xl[j,2]
        mID[ind_shift+j] = mol[j]
        typ[ind_shift+j] = tp[j]
        ind[ind_shift+j] = il[j]



for i in range(nP):
    xl, bl, tp = make_gnp(Rp, nG, Nga, Ngb, box, Dim)
    
    ind_shift = nD*N + nH*Nh + i * ns_per_gnp
    print("ind_shift: " , ind_shift)


    bl[:,1] = bl[:,1] + ind_shift
    bl[:,2] = bl[:,2] + ind_shift
    
    for j in range(bonds_per_gnp):
        bonds[bind,:] = bl[j,:]
        bind = bind + 1


    mol = np.zeros(ns_per_gnp,'i')
    mol[:] = nD+i+1
    
    il = np.zeros(ns_per_gnp,'i')
    for j in range(ns_per_gnp):
        il[j] = ind_shift + j + 1
        x[ind_shift+j,0] = xl[j,0]
        x[ind_shift+j,1] = xl[j,1]
        x[ind_shift+j,2] = xl[j,2]
        mID[ind_shift+j] = i+1 + nD + nH
        typ[ind_shift+j] = tp[j]
        ind[ind_shift+j] = il[j]

 

"""

