/* Lucretia data, prototypes, functions related to mathematics and physics */

/* AUTH: PT, 02-aug-2004 */
/* MOD:
                         */
#ifndef LUCRETIA_COMMON
  #include "LucretiaCommon.h"
#endif
#define LUCRETIA_PHYSICS

/* define some constants */

/* value of pi from Matlab 6.1: */
#define PI 3.14159265358979
/* value of electron mass in GeV from 2004 PDG: */
#define ME2C2_GEV 0.000510998918
/* conversion from E[GeV] to Brho [T.m] */
#define GEV2TM 3.335640952
/* speed of light in m/s */
#define CLIGHT 299792458.0 
/* define energy gain which switches between Lcav matrix and drift mat*/
#define MIN_EGAIN 1.e-12
/* min bunch length, replaces sigz if sigz==0 in SRWF convolution */
#define MIN_BUNCH_LENGTH 0.3e-06
/* 300 nm == 1 femtosecond */
#define MIN_TILT 1.e-12
/* if magnet is within MIN_TILT of being a normal magnet, neglect skew terms */

/* define a macro to compute the fall-behind due to relativistic
   effects (actually the fall-behind is this factor times the element
	length) */

#define LORENTZ_DELAY(P) (ME2C2_GEV *ME2C2_GEV / 2 / P / P)

/* return the LucretiaPhysics version string */ 

char* LucretiaPhysicsVersion( ) ;

/* return the R matrix of a drift space */

void GetDriftMap( double, Rmat ) ;

/* return the R matrix and T5xx terms of a quad */

void GetQuadMap( double, double, double, double, Rmat, double[] ) ;

/* return the R matrix and T5xx terms of a solenoid */

void GetSolenoidMap( double, double, double, Rmat, double[] ) ;

/* return the T matrix terms for a sextupole */

void GetSextMap( double, double, double, double, 
			 double[4][10] ) ;

/* return the R matrix for an RF structure */

void GetLcavMap( double, double, double, double, double,
				 Rmat, int ) ;

/* propagate the transverse coordinates of a ray thru a thin-lens
   multipole */

void PropagateRayThruMult( double, double*, double*, double*, int, double*,
						   double, double, double*, double*, int, int, double,
							double*, int*, double*, double*, int, int ) ;

/* emulation of the MAD transport map for a sector bend */

void GetMADSBendMap( double, double, double, double,
						   double, Rmat, double[4][10], 
							double[10], int ) ;

/* transfer map for a sector bend, Lucretia native form: */

void GetLucretiaSBendMap( double , double, double, double,
						  double, Rmat, double[10] ) ;

/* return the R matrix for a sector bend fringe field */

void GetBendFringeMap( double, double, double, 
					   double, double, double, 
					   double, double, double, 
					   Rmat, double[10] ) ;

/* perform a rotation of an R-matrix through an xy angle */

void RotateRmat( Rmat, double ) ;

/* propagate a set of twiss parameters thru an r matrix */

int TwissThruRmat( Rmat, struct beta0*, struct beta0* ) ;

/* propagate coupled twiss parameters through an R matrix */

int CoupledTwissThruRmat( Rmat, double*, double*, int ) ;

/* convolve a short-range wakefield with the beam */

int ConvolveSRWFWithBeam( struct Bunch*, int, int ) ;

/* bin the rays in a bunch according to the desired bin width */

int BinRays( struct Bunch*, struct SRWF*, double ) ;

/* prepare a bunch to participate in frequency-domain LRWFs */

int PrepareBunchForLRWFFreq(struct Bunch*, int, int,
									 double*, double*, int,
									 double, double ) ;

/* compute the kick from a frequency-domain LRWF */

int ComputeTLRFreqKicks(struct LRWFFreq*, double, int, 
								struct LRWFFreqKick*, int, int,
								struct LucretiaComplex*, 
								struct LucretiaComplex*, 
								struct LucretiaComplex*, 
								double*, double ) ;

/* complex product operation */

struct LucretiaComplex ComplexProduct( struct LucretiaComplex, 
												   struct LucretiaComplex ) ;

/* synchrotron radiation parameters */

void CalculateSRPars( double, double, double, 
							 double*, double*, double*, double* ) ;

/* Poisson-distributed random numbers */

int poidev( double ) ;

/* SR photon distribution via Wolski's method */

double SRSpectrumAW( ) ;

/* SR photon distribution via Burkhardt's method */

double SRSpectrumHB( ) ;
double SynRadC( double ) ;

/* master SR loss function */

double ComputeSRMomentumLoss( double, double, double, int ) ;

/* transfer map for a coordinate-change element */

int GetCoordMap( double[6], double[6], Rmat ) ;
