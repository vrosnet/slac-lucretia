//extern "C" {
  #include "../LucretiaGlobalAccess.h"
  #include "../LucretiaMatlab.h"
//}
#include <string>
#define LUCRETIA_MANAGER

using namespace std;

class lucretiaManager
{
  public:
  lucretiaManager(int* blele, struct Bunch* ThisBunch) ;
  ~lucretiaManager();
  struct Bunch* fBunch ;
  int GetNextX() ;
  void SetNextX(double x[6], int id, int doresume) ;
  int Status ; // 0=OK; 1=No external processes defined for this BEAMLINE element
  char* GeomType ; // Geometry type ("Ellipse" | "Rectangle")
  char* Material ; // Material type -> must be from GEANT4 material tables
  double AperX ; // Horizontal half-aperture / m
  double AperY ; // Vertical half-aperture / m
  double Ecut ; // Energy cut for storing tracks
  double Thickness ; // Material thickness / m
  int Verbose; // GEANT4 text output verbosity
  private:
  int* fEle ;
  int fRayCount ;
  int fRayGetPtr ;
  int* fPrimIndex ;
  double fMaxPrimaryParticles ; // Max number of primaries to generate from Lucretia beam
} ;
