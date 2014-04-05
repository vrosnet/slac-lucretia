#include "actionInitialization.hh"
#include "primaryGeneratorAction.hh"
#include "runAction.hh"
#include "eventAction.hh"
#include "trackingAction.hh"
#ifndef LUCRETIA_MANAGER
  #include "lucretiaManager.hh"
#endif

actionInitialization::actionInitialization(G4double zpos, G4double ecut, G4double zcut, lucretiaManager* lman)
  : G4VUserActionInitialization(),
    fZPOS(zpos),
    fECUT(ecut),
    fZCUT(zcut),
    fLman(lman)
{
}

actionInitialization::~actionInitialization()
{}


void actionInitialization::BuildForMaster() const
{
  SetUserAction(new runAction);
}

void actionInitialization::Build() const
{
  primaryGeneratorAction* prim = new primaryGeneratorAction(fZPOS,fLman) ;
  SetUserAction(prim);
  SetUserAction(new runAction);
  SetUserAction(new eventAction);
  SetUserAction(new trackingAction(fECUT,fZCUT,fLman));
}  

