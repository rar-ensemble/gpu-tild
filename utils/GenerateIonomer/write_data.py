#!/home/rrig/Install/anaconda3/bin/python3
import numpy as np
import os, sys
from numpy.random import rand, normal

def write_data( nm , index , mID , typ, charges, box , x , bonds , angles ):
  nangs = np.shape(angles)[0]
  nbonds = np.shape(bonds)[0]
  ns = np.shape(x)[0]

  do_charges = 0
  for i in range(ns):
      if ( charges[i] != 0.0 ):
          do_charges = 1

  ntypes = np.shape( np.unique(typ) )[0]
  nBondtypes = np.shape( np.unique(bonds[:,0]) )[0]

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

  ln = '%d atom types\n%d bond types\n%d angle types\n\n' % (ntypes,nBondtypes,nang_types)
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
    if ( do_charges ):
        ln = '%d %d %d  %1.3f  %f %f %f\n' % ( index[i] , mID[i] , typ[i] , charges[i], x[i][0] , x[i][1] , x[i][2] )
    else:
        ln = '%d %d %d  %f %f %f\n' % ( index[i] , mID[i] , typ[i] , x[i][0] , x[i][1] , x[i][2] )
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


