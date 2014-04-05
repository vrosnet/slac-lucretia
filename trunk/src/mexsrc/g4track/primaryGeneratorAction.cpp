#include "primaryGeneratorAction.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "math.h"
#ifndef LUCRETIA_MANAGER
  #include "lucretiaManager.hh"
#endif
#include <iostream>

primaryGeneratorAction::primaryGeneratorAction(G4double zpos, lucretiaManager* lman)
  : G4VUserPrimaryGeneratorAction(),
    fZPOS(zpos),  // starting Z postion (left edge of world)
    fLman(lman)
{
  // Input particles type of electron
  G4int nofParticles = 1;
    fParticleGun = new G4ParticleGun(nofParticles);
  G4ParticleDefinition* particleDefinition 
    = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  fParticleGun->SetParticleDefinition(particleDefinition);
  
}

primaryGeneratorAction::~primaryGeneratorAction()
{
  delete fParticleGun;
}

void primaryGeneratorAction::GeneratePrimaries(G4Event* Event)
{
  // This function is called at the begining of event
  // Generate primary particles from Lucretia input (stopped) beam structure
  int xPtr ;
  double Px, Py, Pz, P ;
  //cout << "+-+-+-+-+-+-+- primaryGeneratorAction +-+-+-+-+-+-+-\n" ;
  xPtr = fLman->GetNextX() ;
  while ( xPtr>=0 ) { // Loop through all stopped Lucretia macro-particles
    // X/Y co-ordinates from Lucretia Bunch, place Z at start of geometry
    fParticleGun->SetParticlePosition(G4ThreeVector(fLman->fBunch->x[xPtr*6]*m, fLman->fBunch->x[xPtr*6+2]*m, -fZPOS*m));
    P=fLman->fBunch->x[xPtr*6+5];
    Px=P*sin(fLman->fBunch->x[xPtr*6+1]);
    Py=P*sin(fLman->fBunch->x[xPtr*6+3]);
    Pz=sqrt(P*P-(Px*Px+Py*Py));
    // Momentum directions and energy from Lucretia Bunch
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(Px*GeV,Py*GeV,Pz*GeV));
    fParticleGun->SetParticleEnergy(P*GeV);
    // Create a GEANT4 primary vertex from this macro-particle
    fParticleGun->GeneratePrimaryVertex(Event);
    // Next stopped Lucretia macro-particle
    xPtr = fLman->GetNextX() ;
    //cout << "X/Y/Z : " << fLman->fBunch->x[xPtr*6]*m << " / " << fLman->fBunch->x[xPtr*6+2]*m << " / " << -fZPOS*m << " Px/Py/Pz : " << Px << " / " << Py << " / " << Pz << " P : " << P << "\n" ;
  }
}

