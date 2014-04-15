#ifndef LUCRETIA_MANAGER
  #include "../LucretiaGlobalAccess.h"
#endif
#include "../LucretiaMatlab.h"
#include <string>
#define LUCRETIA_MANAGER
#include "mex.h"

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
  int Status ; // 0=OK; 1=No external processes defined for this BEAMLINE element
  char* GeomType ; // Geometry type ("Ellipse" | "Rectangle")
  char* Material ; // Material type -> must be from GEANT4 material tables
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

  private:
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

