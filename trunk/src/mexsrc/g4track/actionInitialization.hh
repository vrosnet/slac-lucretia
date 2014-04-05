#include "G4VUserActionInitialization.hh"
#include "G4UnitsTable.hh"
#ifndef LUCRETIA_MANAGER
  #include "lucretiaManager.hh"
#endif

//class B4DetectorConstruction;

/// Action initialization class.
///

//class lucretiaManager;

class actionInitialization : public G4VUserActionInitialization
{
public:
  actionInitialization(G4double zpos, G4double ecut, G4double zcut, lucretiaManager* lman );
  ~actionInitialization();
  virtual void BuildForMaster() const;
  virtual void Build() const;
private:
  G4double fZPOS;
  G4double fECUT;
  G4double fZCUT;
  lucretiaManager* fLman ;
};


    
