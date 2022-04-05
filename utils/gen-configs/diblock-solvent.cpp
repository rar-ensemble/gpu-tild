#include <complex>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

using namespace std ;

#define PI   3.141592653589793238462643383
double ran2( void ) ;
double gasdev2( void ) ;
long int idum ;

int main(int argc, char** argv ) {

  idum = int( time(0) ) ;

  int i, m, j, k, ind = 0, n ;

  double L[3], V, C, Rg, rho0, phi;
  int Dim = 2, Nda = 10, Ndb = 0;
  int ncomps = 2 ;
  int *nmolecs = new int [ncomps] ;
  
  L[0] = L[1] = 45.0 ;
  L[2] = 45.0 ;
  
  int *N = new int [ncomps] ;
  N[0] = Nda + Ndb ;
  N[1] = 1;
  
  phi = 0.6;

  // Define C as though this were a polymer melt
  // to allow mapping onto a melt molecular weight
  C = 1.0 ;

  Rg = sqrt( ( double(N[0]) - 1.0 ) / 6.0 ) ;
  V = 1.0 ;
  double RgD = 1.0 ;
  for ( i=0 ; i<Dim ; i++ ) {
    V *= L[i] ;
    RgD *= Rg ;
  }

  
  nmolecs[0] = int( C * phi * V / RgD ) ;
  nmolecs[1] = int( (1.0-phi) * float(nmolecs[0]*N[0]) / phi );
//  cout << "Should agree: " << 1.0-phi << " " << float(nmolecs[1])/float(nmolecs[0]*N[0]+nmolecs[1]) << endl;

  int nstot = 0 ;
  for ( i=0 ; i<ncomps ; i++ )
    nstot += N[i] * nmolecs[i] ;

  double **x = new double*[nstot] ;
  for ( i=0 ; i<nstot; i++ )
    x[i] = new double[Dim] ;

  int *tp = new int [nstot] ;
  int *mID = new int [nstot] ;


  ind = 0 ;
  for ( k=0 ; k<nmolecs[0] ; k++ ) { 
    for ( j=0 ; j<Dim ; j++ )
      x[ind][j] = ran2() * L[j] - L[j]/2.0;
    


    tp[ind] = 0;
    mID[ind] = k ;

    ind++ ;
      
    for ( m=1 ; m<Nda ; m++ ) {
    
      for ( j=0 ; j<Dim ; j++ ) {
        x[ind][j] = x[ ind-1 ][ j ] + gasdev2() ;
 
        if ( x[ind][j] > L[j]/2.0 )
          x[ind][j] -= L[j] ;
        else if ( x[ind][j] < -L[j]/2.0 )
          x[ind][j] += L[j] ;
      }
 
      tp[ ind ] = 0;
      mID[ind] = k ;
 
      ind++ ;
      
    }
    for ( m=Nda ; m<Nda+Ndb ; m++ ) {
    
      for ( j=0 ; j<Dim ; j++ ) {
        x[ind][j] = x[ ind-1 ][ j ] + gasdev2() ;
 
        if ( x[ind][j] > L[j]/2.0 )
          x[ind][j] -= L[j] ;
        else if ( x[ind][j] < -L[j]/2.0 )
          x[ind][j] += L[j] ;
      }
 
      tp[ ind ] = 1;
      mID[ind] = k ;
 
      ind++ ;
      
    }
  }// for n=0:nmolecs[0];

  // PLACE THE SOLVENT AT RANDOM //
  for ( k=nmolecs[0] ; k<nmolecs[0]+nmolecs[1] ; k++ ) { 
    for ( j=0 ; j<Dim ; j++ )
      x[ind][j] = ran2() * L[j] - L[j]/2.0;

    tp[ind] = 2;
    mID[ind] = k ;

    ind++ ;
  }


  rho0 = double(nmolecs[0]*N[0]+nmolecs[1]) / V ;
  C = nmolecs[0] / V * RgD ;
  ofstream otp("input.data") ;

  otp << "C: " << C << " rho0: " << rho0 << " [b^-Dim], phiS: " <<
    phi << "\n\n" ;
  otp << nstot << " atoms\n" ;
  otp << nmolecs[0]*(N[0]-1) << " bonds\n";
  otp << "0 angles\n\n" ;

  otp << "3 atom types\n1 bond types\n0 angle types\n\n" ; 

  otp << -L[0]/2.0 << " " << L[0]/2.0 << " xlo xhi\n" ;
  otp << -L[1]/2.0 << " " << L[1]/2.0 << " ylo yhi\n" ;
  if ( Dim > 2 )
    otp << -L[2]/2.0 << " " << L[2]/2.0 << " zlo zhi\n\n" ;
  else
    otp << "-5.0 5.0 zlo zhi\n\n" ;

  otp << "Masses\n\n" ; 
  otp << "1 1.0\n2 1.0\n3 1.0\n\n" ;

  otp << "Atoms\n\n" ;

  for ( i=0 ; i<nstot; i++ ) {
    otp << i+1 << " " << mID[i]+1 << " " << tp[i]+1 << " " ;

    for ( j=0 ; j<Dim ; j++ ) 
      otp << x[i][j] << " " ;
    
    for ( j=Dim ; j<3 ; j++ )
      otp << "0.0" ;
    
    otp << endl;
  }

  otp << "\nBonds\n\n" ;
  
  int bind = 0 ;
  ind = 1 ;
  for ( i=0 ; i<ncomps ; i++ ) {
    for ( j=0 ; j<nmolecs[i] ; j++ ) {
      for ( k=0 ; k<N[i]-1 ; k++ ) {
        otp << ind << " 1 " << j*N[i]+k+1 << " " << j*N[i]+k+2 << endl;
        ind++ ;
      }
    }
  }

  otp.close() ;

  printf("Random config generated!\n") ; fflush( stdout ) ;
  return 0 ;

}

/* ======================================================================== */
/* random.c                                                                 */
/*            Random number generator from uniform distribution [0,1]       */
/*             from Numerical Recipes                                       */
/* ======================================================================== */
#include <cmath>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0 / IM1)
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1 + IMM1 / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)


double ran2 (void)
{
  int j;
  long int k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  double temp;
  extern long idum;
  
  if (idum <= 0) {
    if (-(idum) < 1) idum = 1;
    else idum = -(idum);
    idum2 = (idum);
    for (j = NTAB + 7; j >=0; j--) {
      k = (idum) / IQ1;
      idum = IA1 * (idum - k * IQ1) - k * IR1;
      if (idum < 0) idum += IM1;
      if (j < NTAB) iv[j] = idum;
    }
    iy = iv[0];
  }
  
  k = (idum) / IQ1;
  idum = IA1 * (idum - k * IQ1) - k * IR1;
  if (idum < 0) idum += IM1;
  k = idum2 / IQ2;
  idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
  if (idum2 < 0) idum2 += IM2;
  j = iy / NDIV;
  iy = iv[j] - idum2;
  iv[j] = idum;
  if (iy < 1) iy += IMM1;
  if ((temp = AM * iy) > RNMX) return RNMX;
  else return temp;
}


void random_unit_vec( double *v ) {
  double xi[3], mxi, mxi2;
  int i;

  mxi2 = 23432.;

  while ( mxi2 > 1.0 ) {
    mxi2 = 0.0 ;
    for ( i=0 ; i<3 ; i++ ) {
      xi[i] = 1.0 - 2.0 * ran2() ;
      mxi2 += xi[i] * xi[i] ;
    }

  }

  mxi = sqrt( mxi2 );

  for ( i=0 ; i<3 ; i++ )
    v[i] = xi[i] / mxi ;
} 


double gasdev2() {
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;
  if(iset == 0) {
    do {
      v1=2.0*ran2()-1.0;
      v2=2.0*ran2()-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);

    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;

  } 
  else {
    iset=0;
    return gset;
  }
}

