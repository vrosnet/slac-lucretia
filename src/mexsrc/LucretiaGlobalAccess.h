/* LucretiaGlobalAccess.h:
   This is where data definitions and function prototypes are put which
   relate to accessing the Lucretia global data objects (notably the
   BEAMLINE, GIRDER, KLYSTRON, PS, WAKEFIELD objects).  The global access
   routines live in LucretiaMatlab.c or LucretiaOctave.c.  However, they
   need to be accessed by procedures in LucretiaCommon.c.  Therefore the
   information is stored here and not in LucretiaCommon.h, nor in either
   LucretiaMatlab.h nor LucretiaOctave.h. */

/* AUTH:  PT, 03-aug-2004. */
/* MOD:                  
        08-mar-2006, PT:
		     add GetMatrixNormalizer prototype.
                           */

#define LUCRETIA_GLOBAL_ACCESS

/* define an enumeration type for klystron status */

enum KlystronStatus {ON, STANDBY, TRIPPED, MAKEUP, STANDBYTRIP} ;

/* function which returns the element class ("DRIF","QUAD",etc) */

char* GetElemClass( int ) ;

/* function to get a numeric parameter from a beamline element: */

double* GetElemNumericPar( int, char*, int* ) ;

/* function to get the index of the power supply for a magnet */

int GetPS( int ) ;

/* get a numeric parameter from a power supply */

double* GetPSNumericPar( int, char*, int* ) ;

/* function to get the index of the girder for an element */

int GetElemGirder( int ) ;

/* Get a numeric paramter for a girder */

double* GetGirderNumericPar( int, char*, int* ) ;

/* get the klystron number for an element */

int GetKlystron( int ) ;

/* get a klystron numeric parameter */

double* GetKlystronNumericPar( int, char*, int* ) ;

/* get a klystron's status */

enum KlystronStatus* GetKlystronStatus( int ) ;

/* how many wakefields of each type */

int* GetNumWakes( ) ;

/* tells how many elements are defined */

int nElemInBeamline( ) ;

/* get other global data sizes */

int GetnGirder( ) ;
int GetnKlystron( ) ;
int GetnPS( ) ;


/* return the total number of track flags an element has */

int GetNumTrackFlags( int elemno ) ;

/* return the name and value of one tracking flag */

int GetTrackFlagValue( int ) ;
const char* GetTrackFlagName( int ) ;

/* Add an error message to the pile */

void AddMessage( const char*, int ) ;

/* retrieve the messages and clear the pile */

char** GetAndClearMessages( int* ) ;

/* use Matlab randn function to get a vector of Gaussian-
   distributed random numbers */

double* RanGaussVecPtr( int ) ;

/* use Matlab rand function to get a vector of uniform-
   distributed random numbers */

double* RanFlatVecPtr( int ) ;

/* Use Matlab sort function to get a sortkey for the rays in  
   a bunch, along a given DOF */

double* GetRaySortkey( double*, int, int ) ;

/* Access the WF global and return pointers to the z positions, kick
   factors; return the bin width as well */

int GetSRWFParameters( int, int, double**, double**, double* ) ;

/* cubic-spline the SRWF and return the splined values at requested
   z locations */

double* SplineSRWF( double*, double*, double*, int, int, int ) ;

/* Use Matlab's pascal function to get a Pascal matrix for use in
   multipole field expansions */

double* GetPascalMatrix( ) ;
double* GetFactorial( ) ;
double  GetMaxMultipoleIndex( ) ;
void    ClearMaxMultipoleStuff( ) ;
void    ComputeNewMultipoleStuff( double ) ;

/* Get the class of a transverse-long-range wakefield (ie, time or
   frequency domain */

int GetTLRWakeClass( int ) ;

/* Get the class of a transverse-long-range error wakefield (ie, time or
   frequency domain */

int GetTLRErrWakeClass( int ) ;

/* get a numeric parameter from a TLR */

double* GetTLRNumericPar( int, char*, int* ) ;

/* get a numeric parameter from an error TLR */

double* GetTLRErrNumericPar( int, char*, int* ) ;

/* get the shape parameter for a collimator */

int GetCollimatorGeometry( int ) ;

/* get log of gamma function */

double GammaLog( double ) ;

/* get cube root of R-matrix determinant */

double GetMatrixNormalizer( double* ) ;

/* Calculation and application of CSR wake */

void GetCsrEloss(struct Bunch*, int, int, int, double, double ) ;


