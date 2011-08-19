/* Data types and function prototypes which are used in the Matlab
   implementation of Lucretia but not in the Octaver version. */

/* AUTH: PT, 02-aug-2004 */
/* MOD:
		   06-Mar-2008, PT:
			   add IsEmpty proto.
		   08-mar-2006, PT:
			   add GetTwissCoupledSetReturn prototype.
                         */

#ifndef mex_h
  #include "mex.h"
#endif
#ifndef LUCRETIA_COMMON
  #include "LucretiaCommon.h"
#endif

/* for GetRmats:  check arguments and set returns */

struct RmatArgStruc* GetRmatsGetCheckArgs( int, int, const mxArray*[] ) ;
void   GetRmatsSetReturn( int, mxArray*[], const struct Rstruc* ) ;

/* for RmatAtoB:  check arguments and set returns */

struct RmatArgStruc* RmatAtoBGetCheckArgs( int, int, const mxArray*[] ) ;
void   RmatAtoBSetReturn( int, mxArray*[], double* ) ;

/* for GetTwiss:  check arguments, unpack twiss, set returns */

struct RmatArgStruc* GetTwissGetCheckArgs( int, int, const mxArray*[] ) ;
double *GetTwissUnpackTwiss( const mxArray* ) ; 
void   ReverseTwissDirection( struct RmatArgStruc* ) ;
void   GetTwissSetReturn( int, mxArray*[], struct twiss* ) ;
void   GetTwissCoupledSetReturn( int, int, mxArray*[], 
										  struct Ctwiss * ) ;

/* for TrackThru:  check arguments and set returns */

struct TrackArgsStruc* TrackThruGetCheckArgs( int, mxArray*[],
														    int, const mxArray*[] ) ;
void TrackThruSetReturn( struct TrackArgsStruc*, int, mxArray*[], int ) ;

/* version information */

mxArray* LucretiaMatlabVersions( const char* ) ;

/* setup global data structure access */

int LucretiaMatlabSetup( ) ;

/* copy real values into a field of a structure */

int SetDoublesToField( mxArray*, int, char*, double*, int ) ;

/* construct the cell matrix for the status + messages */

mxArray* CreateStatusCellArray( int, int, char** ) ;

/* general purpose tool for divining the class of a long-range wake */

int GetTLRGeneralWakeClass( int, const mxArray* ) ;

/* general purpose tool for interrogating the lrwf database */

double* GetTLRGeneralNumericPar( int, char*, int*, const mxArray* ) ;

/* clear return variables from randn and spline calls */

void ClearPLHSVars( ) ;

/* Check in a version-neutral way to see whether an mxArray is empty */

bool IsEmpty( const mxArray* ) ;

/* Crazily enough, Mathworks uses mxCreateDoubleScalar for version 6.5 and
   later, but mxCreateScalarDouble for 6.1 and earlier.  There is a similar
   issue with mexGetArrayPtr and mexGetVariablePtr, except that here there
   is a change in the variable order!
   
	 In order to 
	finesse this, we offer the following redefinitions: */

/* for earlier versions */

/* #define CREATE_SCALAR_DOUBLE mxCreateScalarDouble
   #define GET_GLOBAL_PTR(A,B) mexGetArrayPtr( A, B ) */

/* for later versions */

#define CREATE_SCALAR_DOUBLE mxCreateDoubleScalar 
#define GET_GLOBAL_PTR(A,B) mexGetVariablePtr( B, A )


