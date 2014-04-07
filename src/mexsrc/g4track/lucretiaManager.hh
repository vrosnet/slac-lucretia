#ifndef LUCRETIA_MANAGER
  #include "../LucretiaGlobalAccess.h"
#endif
#include "../LucretiaMatlab.h"
#include <string>
#define LUCRETIA_MANAGER

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
  void SetNextSecondary(double x[6], const char* type) ;
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
  double fMaxSecondaryParticles ; // Max number of secondaries to keep
  double fMaxSecondaryParticlesPerPrimary ; // Max secondaries to store per primary launched
  int fSecondariesCounter ;
  int* fSecondariesPerThisPrimary ;

  private:
  double* fNumSecondariesStored ;
  double* fSecondaryBunch_x ;
  double* fSecondaryBunch_q ;
  mxArray* fSecondaryBunch_type ;
  int* fBunchNo ; // Bunch number
  double* fPrimOrder ; // Order in which to read in primary particles
  int* fEle ;
  int fRayCount ;
  int fRayGetPtr ;
  int* fPrimIndex ;
  double fMaxPrimaryParticles ; // Max number of primaries to generate from Lucretia beam
} ;

