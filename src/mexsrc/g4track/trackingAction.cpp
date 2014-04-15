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

trackingAction::trackingAction(lucretiaManager* lman)
 :G4UserTrackingAction(),
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
  G4double momx = track->GetMomentum().x()/CLHEP::GeV ;
  G4double momy = track->GetMomentum().y()/CLHEP::GeV ;
  G4double momz = track->GetMomentum().z()/CLHEP::GeV ;
  G4int parentID = track->GetParentID() ;
  double xLucretia[6] ;
  int passCuts=0 ;
  int dosecondaries = 0;
  static int hitsPerParentCounter ;
  static int lastPrimaryID ;
  if ( fLman->fMaxSecondaryParticles>0 && fLman->fMaxSecondaryParticlesPerPrimary>0)
    dosecondaries = 1 ;
  if (parentID==0 || dosecondaries ) { 
    xLucretia[0] = x; xLucretia[2] = y; xLucretia[4] = z; // z doesn't get copied back to Lucretia bunch
    xLucretia[1] = atan(momx/momz) ; xLucretia[3] = atan(momy/momz) ;
    xLucretia[5] = e ;
    if (z>=fLman->Lcut && momz>=fLman->Ecut && e>=fLman->Ecut)
      passCuts=1;
    if (dosecondaries && passCuts &&
            fLman->fSecondariesCounter < fLman->fMaxSecondaryParticles )
      dosecondaries = 1 ;
    else
      dosecondaries = 0 ;
  }
  if (parentID==0) { // Actions for primary particles
    hitsPerParentCounter=0;
    //std::cout << "Z: " << z << " Z_cut: " << fLman->Lcut << " momz: " << momz << " E: " << e << " E_cut: " << fLman->Ecut << " CUTS: " << passCuts << "\n" ;
    fLman->SetNextX(xLucretia, track->GetTrackID()-1, passCuts) ; // Write new tracked particle back to Lucretia Matlab Bunch
    lastPrimaryID=track->GetTrackID()-1;
    /*std::cout << "-=-=-=-=-***USER_TRACKING****-=-=-=-= [PRIMARY PARTICLE]\n"
	      << "ID= " << track->GetTrackID() << "  e= " << e << " Name:  "
	      << track->GetDynamicParticle()->GetDefinition()->GetParticleName() << "  "
	      << "X/Y/Z: " << x << " / " << y << " / " << z << "  "
	      << "PX/PY: " << atan(momx/momz) << " / " << atan(momy/momz) << "  "
	      << "P: " << momx << " , " << momy << " , " << momz << " "
	      << "Parent ID= " << parentID
	      << "\n"  ;*/
  }
  else if (dosecondaries && (uint32_T)hitsPerParentCounter<fLman->fMaxSecondaryParticlesPerPrimary) {
    hitsPerParentCounter++;
    /*std::cout << "-=-=-=-=-***USER_TRACKING****-=-=-=-= [SECONDARY PARTICLE]\n"
	      << "ID= " << track->GetTrackID() << "  e= " << e << " Name:  "
	      << track->GetDynamicParticle()->GetDefinition()->GetParticleName() << "  "
	      << "X/Y/Z: " << x << " / " << y << " / " << z << "  "
	      << "PX/PY: " << atan(momx/momz) << " / " << atan(momy/momz) << "  "
	      << "P: " << momx << " , " << momy << " , " << momz << " "
	      << "Parent ID= " << parentID
	      << "\n"  ;*/
    fLman->SetNextSecondary(xLucretia, lastPrimaryID, track->GetDynamicParticle()->GetDefinition()->GetParticleName()) ; 
  }
}
