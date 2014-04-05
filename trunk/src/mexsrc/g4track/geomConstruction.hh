#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "CLHEP/Units/SystemOfUnits.h"

class G4LogicalVolume;
class G4Material;
class DetectorMessenger;

class geomConstruction : public G4VUserDetectorConstruction
{
  public:
    geomConstruction(
      const G4String& collType = "Rectangle",
      const G4String& collMaterialName = "G4_Be",
      G4double aper_x = 0*CLHEP::m,
      G4double aper_y = 0*CLHEP::m,
      G4double thickness = 1*CLHEP::m,
      G4double dz = 1*CLHEP::m);
    ~geomConstruction();

  public:
    // methods from base class 
    virtual G4VPhysicalVolume* Construct();

  private:
  G4String fCollType;
  G4String fCollMaterialName;
  G4double fCollAperX;
  G4double fCollAperY;
  G4double fCollThickness;
  G4double fCollLength;
};

