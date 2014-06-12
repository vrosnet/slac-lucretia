#ifndef LUCRETIA_MANAGER
  #include "../LucretiaGlobalAccess.h"
#endif
#include "../LucretiaMatlab.h"
#include <string>
#define LUCRETIA_MANAGER
#include "mex.h"

/* Interpolation code extracted from modified ba_interp3.cpp from Mathworks File Exchange*/
/* ba_inter3 GPL License (c) 2008 Brian Amberg http://www.brian-amberg.de/ */

using namespace std;

class lucretiaManager
{
  public:
  lucretiaManager(int* blele, int* bunchno, struct Bunch* ThisBunch, double L) ;
  ~lucretiaManager();
  void Initialize(int* blele, int* bunchno, struct Bunch* ThisBunch, double L) ;
  struct Bunch* fBunch ;
  int GetNextX() ;
  void SetNextX(double x[6], int id, int doresume) ;
  void SetNextSecondary(double x[6], int id, const char* type) ;
  void SetLucretiaData() ;
  void freeMem() ;
  void GetUniformField(double uField[3]) ;
  double interpField(const int fieldno, const double* point) ;
  int Status ; // 0=OK; 1=No external processes defined for this BEAMLINE element
  char* GeomType ; // Geometry type ("Ellipse" | "Rectangle" | "Tapered")
  char* Material ; 
  char* Material2 ; 
  char* VacuumMaterial ;
  const mxArray *pBx,*pBy,*pBz,*pEx,*pEy,*pEz ; // EM field values
  char* EMStepperMethod ;
  double EMStepSize ;
  double EMDeltaOneStep ;
  double EMDeltaIntersection ;
  double EMDeltaChord ;
  double EMEpsMin ;
  double EMEpsMax ;
  char* EMInterpMethod ;
  char EnableEM ;
  char EMisUniform ;
  double AperX ; // Horizontal half-aperture / m
  double AperY ; // Vertical half-aperture / m
  double Ecut ; // Energy cut for storing tracks
  double Thickness ; // Material thickness / m
  int Verbose; // GEANT4 text output verbosity
  double Lcut; // Final track must be >= this to make it back into Lucretia
  int fNumRaysResumed ; // Keep count of number of un-stopped rays
  uint32_T fMaxSecondaryParticles ; // Max number of secondaries to keep
  uint32_T fMaxSecondaryParticlesPerPrimary ; // Max secondaries to store per primary launched
  uint32_T fSecondariesCounter ;
  int* fSecondariesPerThisPrimary ;
  int* fEle ;
  double CollLen2 ;
  double AperX2 ;
  double AperY2 ;
  double AperX3 ;
  double AperY3 ;
  double CollDX ;
  double CollDY ;
  struct UserElement {
    char* Name;
    char* Symbol;
    double Z;
    double A;
    double FractionMass;
  } ;
  struct {
	  double density;
    double pressure;
    double temperature;
    char* state;
    size_t NumComponents;
    struct UserElement* element;
  } UserMaterial[3];

  private:
  int access(int M, int N, int O, int x, int y, int z) ;
  int access_unchecked(int M, int N, int O, int x, int y, int z) ;
  void indices_linear(
        int &f000_i,
        int &f100_i,
        int &f010_i,
        int &f110_i,
        int &f001_i,
        int &f101_i,
        int &f011_i,
        int &f111_i,
        const int x, const int y, const int z,
        const mwSize &M, const mwSize &N, const mwSize &O) ;
  void indices_cubic(int f_i[64], const int x, const int y, const int z,
                                      const mwSize &M, const mwSize &N, const mwSize &O) ;
  void interpolate_nearest(double *pO, const double *pF,
        const double *pX, const double *pY, const double *pZ,
        const mwSize ND, const mwSize M, const mwSize N, const mwSize O, const mwSize P,
        const double s_x, const double o_x,
        const double s_y, const double o_y,
        const double s_z, const double o_z) ;
  void interpolate_linear(double *pO, const double *pF,
        const double *pX, const double *pY, const double *pZ,
        const mwSize ND, const mwSize M, const mwSize N, const mwSize O, const mwSize P,
        const double s_x, const double o_x,
        const double s_y, const double o_y,
        const double s_z, const double o_z) ;
  void interpolate_bicubic(double *pO, const double *pF,
        const double *pX, const double *pY, const double *pZ,
        const mwSize ND, const mwSize M, const mwSize N, const mwSize O, const mwSize P,
        const double s_x, const double o_x,
        const double s_y, const double o_y,
        const double s_z, const double o_z) ;
  mxArray* fTypeCellPtr ;
  uint32_T* fPrimaryRegenID ;
  uint32_T* fSecondaryPrimaryID ;
  double* fSecondaryBunch_x ;
  int* fBunchNo ; // Bunch number
  uint32_T* fPrimOrder ; // Order in which to read in primary particles
  uint32_T fRayCount ;
  int fRayGetPtr ;
  int* fPrimIndex ;
  uint32_T fMaxPrimaryParticles ; // Max number of primaries to generate from Lucretia beam
} ;

