#include "trackingAction.hh"
#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
#include <math.h>
#include <iostream>
#include "CLHEP/Units/SystemOfUnits.h"
#ifndef LUCRETIA_MANAGER
  #include "lucretiaManager.hh"
#endif

trackingAction::trackingAction(G4double ecut, G4double zcut, lucretiaManager* lman)
 :G4UserTrackingAction(),
  fECUT(ecut),
  fZCUT(zcut),
  fLman(lman)
{ }


void trackingAction::PreUserTrackingAction(const G4Track*)
{
}

void trackingAction::PostUserTrackingAction(const G4Track* track)
{
  G4ThreeVector pos = track->GetPosition() ;
  G4double x = track->GetPosition().x()/CLHEP::m ;
  G4double y = track->GetPosition().y()/CLHEP::m ;
  G4double z = track->GetPosition().z()/CLHEP::m ;
  G4double e = track->GetDynamicParticle()->GetKineticEnergy()/CLHEP::GeV ;
  G4double momx = track->GetMomentum().x() ;
  G4double momy = track->GetMomentum().y() ;
  G4double momz = track->GetMomentum().z() ;
  G4int parentID = track->GetParentID() ;
  double xLucretia[6] ;
  int resumeTracking=0 ;
  if (parentID==0 ) { // Actions for primary particles
    xLucretia[0] = x; xLucretia[2] = y; xLucretia[4] = z; // z doesn't get copied back to Lucretia bunch
    xLucretia[1] = atan(momx/momz) ; xLucretia[3] = atan(momy/momz) ;
    xLucretia[5] = e ;
    if (z>=fZCUT && momz>=fZCUT && e>=fECUT)
      resumeTracking=1;
    fLman->SetNextX(xLucretia, track->GetTrackID()-1, resumeTracking) ; // Write new tracked particle back to Lucretia Matlab Bunch
    /*std::cout << "-=-=-=-=-***USER_TRACKING****-=-=-=-= [PRIMARY PARTICLE]\n"
	      << "ID= " << track->GetTrackID() << "  e= " << e << " Name:  "
	      << track->GetDynamicParticle()->GetDefinition()->GetParticleName() << "  "
	      << "X/Y/Z: " << x << " / " << y << " / " << z << "  "
	      << "PX/PY: " << atan(momx/momz) << " / " << atan(momy/momz) << "  "
	      << "P: " << momx << " , " << momy << " , " << momz << " "
	      << "Parent ID= " << parentID
	      << "\n"  ;*/
  }
}
