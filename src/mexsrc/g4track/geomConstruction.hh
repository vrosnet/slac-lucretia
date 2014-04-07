#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4EllipticalTube.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh" 
#include "G4Box.hh"
#include "G4NistManager.hh"
#ifndef LUCRETIA_MANAGER
  #include "lucretiaManager.hh"
#endif
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
    // others
  void SetGeomParameters(lucretiaManager* lman);

  private:
  G4String fCollType;
  G4String fCollMaterialName;
  G4double fCollAperX;
  G4double fCollAperY;
  G4double fCollThickness;
  G4double fCollLength;
  G4NistManager* nistManager ;
  G4Element* elN;
  G4Element* elO;
  G4Material* Vacuum;
  G4Box* sWorld;
  G4LogicalVolume* worldVolume;
  G4VPhysicalVolume* pWorld;
  G4Box* collBox;
  G4EllipticalTube* collTube ;
  G4LogicalVolume* collVolume;
  G4PVPlacement* collPlacement;
  G4Box* collInnerBox ;
  G4EllipticalTube* collInnerTube ;
  G4LogicalVolume* collVolumeInner ;
  G4PVPlacement* collInnerPlacement ;
};

