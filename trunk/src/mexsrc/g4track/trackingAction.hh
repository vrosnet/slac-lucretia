#include "G4UserTrackingAction.hh"
#include "G4UnitsTable.hh"
#ifndef LUCRETIA_MANAGER
  #include "lucretiaManager.hh"
#endif

class trackingAction : public G4UserTrackingAction {

  public:
  trackingAction(G4double ecut, G4double zcut, lucretiaManager* lman);
  virtual ~trackingAction() {};

  virtual void  PreUserTrackingAction(const G4Track*);
  virtual void PostUserTrackingAction(const G4Track*);
private:
  G4double fECUT;
  G4double fZCUT;
  lucretiaManager* fLman ;
};
