#!/home/rrig/Install/anaconda3/bin/python3
import numpy as np
import os, sys
from numpy.random import rand, normal

def write_data( nm , index , mID , type , box , x , bonds , angles ):
  nangs = np.shape(angles)[0]
  nbonds = np.shape(bonds)[0]
  ns = np.shape(x)[0]

  ntypes = np.shape( np.unique(type) )[0]

  otp = open( nm , 'w' )

  otp.write( 'arrrr\n\n' )

  ln = '%d atoms\n' % ns
  otp.write( ln )

  ln = '%d bonds\n' % nbonds
  otp.write( ln )

  ln = '%d angles\n\n' % nangs
  otp.write( ln )

  nang_types = 0
  if ( nangs > 0 ):
    nang_types = 1

  ln = '%d atom types\n%d bond types\n%d angle types\n\n' % (ntypes,1,nang_types)
  otp.write( ln )

  ln = '%f %f xlo xhi\n' % (0.0 , box[0])
  otp.write( ln )

  ln = '%f %f ylo yhi\n' % (0.0 , box[1])
  otp.write( ln )

  ln = '%f %f zlo zhi\n\n' % (0.0 , box[2])
  otp.write( ln )

  ln = 'Masses\n\n'
  otp.write( ln )
  for i in range(ntypes):
      ln = '%d %f\n' % (i+1, 1.0)
      otp.write(ln)
  otp.write('\n')

  otp.write( 'Atoms\n\n' )

  for i in range ( 0 , ns ):
    ln = '%d %d %d  %f %f %f\n' % ( index[i] , mID[i] , type[i] , x[i][0] , x[i][1] , x[i][2] )
    otp.write( ln )


  otp.write( '\nBonds\n\n' )
  for i in range( 0 , nbonds ):
    ln = '%d %d  %d %d\n' % (i+1 , bonds[i][0] , bonds[i][1] , bonds[i][2])
    otp.write( ln )

  if ( nangs > 0 ):

    otp.write( 'Angles\n\n' )
    for i in range( 0 , nangs ):
      ln = '%d %d  %d %d %d\n' % ( i , angles[i][0] , angles[i][1] , angles[i][2] , angles[i][3])
      otp.write( ln )


  otp.close()


def make_gnp(Rp, nG, Nga, Ngb, box, Dim):
    N = Nga + Ngb
    ns = nG*(N+1)+1

    r = np.zeros((ns,3),'d')
    bonds = np.zeros((nG*N,3),'i')
    typ = np.zeros(ns,'i')

    for j in range(Dim):
        r[0,j] = rand() * box[j] 
    typ[0] = 3

    index = 1
    bindex = 0
    for i in range(nG):
        rg = np.zeros(3,'d')
        for j in range(Dim):
            rg[j] = normal()
        rg = rg * Rp / np.sum(rg**2)**0.5

        for j in range(Dim):
            r[index,j] = r[0,j] + rg[j]
        
        typ[index] = 4
        
        bonds[bindex,0] = 1
        bonds[bindex,1] = index + 1
        bonds[bindex,2] = index + 2

        index = index + 1
        bindex = bindex + 1


        
        for i in range(N):
            for j in range(Dim):
                r[index,j] = r[0,j] + rg[j]/Rp*(Rp + 0.5*(i+1))
            if ( i < Nga ):
                typ[index] = 1
            else:
                typ[index] = 2

            if ( i < N-1 ):
                bonds[bindex,0] = 1
                bonds[bindex,1] = index + 1
                bonds[bindex,2] = index + 2
                bindex = bindex + 1

            index = index + 1
    return r, bonds, typ


# Discrete nanoparticles
def make_one_dnp(Rp, rho0, box, Dim):
    vp = 4. * np.pi * Rp**2
    if ( Dim == 3 ):
        vp = vp * Rp / 3.0

    approx_n = rho0 * vp 
    #print("n roughly ", approx_n)

    d = (vp/approx_n)**(1.0/Dim)

    #print('approx spacing:', d)

    x = np.array([])
    y = np.array([])
    z = np.array([])

    n_r_points = int(2.0*Rp/d)
    ro = np.linspace(-Rp,Rp,n_r_points)
    #ro = np.arange(-Rp,Rp+d,d)
    rz = np.array([0.0])
    if ( Dim == 3 ):
        rz = ro

    nkept = 0
    for xl in ro:
        for yl in ro:
            for zl in rz:
                magr = xl*xl + yl*yl + zl*zl
                if (magr <= Rp+d):
                    x = np.append(x,xl)
                    y = np.append(y,yl)
                    z = np.append(z,zl)
                    nkept += 1

                
    print('nsites:', nkept, 'nsites/vp:', nkept/vp)

    min_bonds = 23483223
    max_bonds = -328923
    bonds = np.zeros(nkept,'i')
    bondEq = np.array([])

    nbonds = 0
    for i in range(nkept):
        for j in range(nkept):
            mdr = ( (x[i] - x[j])**2 + (y[i]-y[j])**2 + (z[i]-z[j])**2)**0.5
            if ( mdr <= 2.1 * d and mdr > 0.):
                bonds[i] = bonds[i] + 1
                if ( j > i ):
                    nbonds = nbonds + 1
                
                #if ( np.count_nonzero(bondEq == mdr) == 0 ):
                if ( np.count_nonzero(np.fabs(bondEq - mdr) < 0.001) == 0 ):
                    bondEq = np.append(bondEq, mdr)
                
            
        if ( bonds[i] > max_bonds):
            max_bonds = bonds[i]
            
        if ( bonds[i] < min_bonds):
            min_bonds = bonds[i]
            
    blist = np.zeros((nbonds,3),'i')
    b_ind = 0
    for i in range(nkept):
        for j in range(i+1,nkept):
            mdr = ( (x[i] - x[j])**2 + (y[i]-y[j])**2 + (z[i]-z[j])**2)**0.5
            if ( mdr <= 2.1 * d and mdr > 0.):
                b_typ = np.where(np.fabs(bondEq-mdr) < 0.001)[0][0]
                blist[b_ind,0] = b_typ + 2
                blist[b_ind,1] = i + 1
                blist[b_ind,2] = j + 1
                b_ind = b_ind + 1

    r = np.column_stack((xl,yl,zl))
    typ = np.zeros(nkept,'i')
    typ[:] = 3

    return r, blist, typ



def make_diblock(Na, Nb, box, Dim):
    typ = np.zeros(Na+Nb,'i')
    r = np.zeros((Na+Nb,3),'d')
    bonds = np.zeros((Na+Nb-1,3),'i')

    for j in range(0,Dim):
        r[0,j] = rand() * box[j]

    if ( Na > 0 ):
        typ[0] = 1
    else:
        typ[0] = 2

    bonds[0,0] = 1;
    bonds[0,1] = 1;
    bonds[0,2] = 2;

    for i in range(1,Na+Nb):
        for j in range(0,Dim):
            r[i,j] = r[i-1,j] + normal();
            if ( r[i,j] > box[j] ):
                r[i,j] = r[i,j] - box[j]
            elif (r[i,j] < 0.0 ):
                r[i,j] = r[i,j] + box[j]

        if ( i < Na ):
            typ[i] = 1
        else:
            typ[i] = 2

        if ( i < (Na+Nb-1)):
            bonds[i,0] = 1;
            bonds[i,1] = i+1;
            bonds[i,2] = i+2;



    return r, bonds, typ
# End of make_diblock routine



box = np.zeros(3,'d')
Dim = 2
rho0 = 4.476
phiD = 1.0
Na = 3
Nb = 6
box[0] = box[1] = 35.0
box[2] = 5.0 ;#2.0*(Na+Nb)**0.5/6.0**0.5
name = "input.data"


N = Na + Nb
V = 1.0
for j in range(0,Dim):
    V = V * box[j]

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
    xl, bl, tp = make_diblock(Na, Nb, box, Dim)
    
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

 


write_data(name, ind, mID, typ, box, x, bonds, [])
