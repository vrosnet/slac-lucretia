#include "geomConstruction.hh"
#include "G4Material.hh" 
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4EllipticalTube.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>

geomConstruction::geomConstruction(const G4String& collType, const G4String& materialName,
				   G4double aper_x, G4double aper_y, G4double thickness, G4double dz)
: G4VUserDetectorConstruction(),
  fCollType(collType),
  fCollMaterialName(materialName),
  fCollAperX(aper_x),
  fCollAperY(aper_y),
  fCollThickness(thickness),
  fCollLength(dz)
{
}
    
  
geomConstruction::~geomConstruction()
{
}

G4VPhysicalVolume* geomConstruction::Construct()
{
  // Define materials via NIST manager
  //
  G4NistManager* nistManager = G4NistManager::Instance();

  // Make Vacuum material for World volume
  G4double density,a,z;
  G4double temperature, pressure, fractionmass;
  G4String name, symbol;
  G4int ncomponents;
  density = 1.290*mg/cm3;
  pressure    = 1.e-39*pascal;
  temperature = 0.11*kelvin;
  a = 14.01*g/mole;
  G4Element* elN  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);
  a = 16.00*g/mole;
  G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);
  G4Material* Vacuum = new G4Material(name="Vacuum"  , density, ncomponents=2, kStateGas,
					   temperature, pressure);
  Vacuum->AddElement(elN, fractionmass=0.7);
  Vacuum->AddElement(elO, fractionmass=0.3);

  // Collimator material
  G4Material* collMaterial
    = nistManager->FindOrBuildMaterial(fCollMaterialName);

  // Geometry parameters
  //
  G4ThreeVector worldDimensions = G4ThreeVector((fCollAperX+fCollThickness)*m,(fCollAperY+fCollThickness)*m,fCollLength*m);

  //std::cout << "WORLD DIMENSIONS: " << worldDimensions.x()/m << " / " << worldDimensions.y()/m << " / " <<  worldDimensions.z()/m << "\n" ;
  
  // World
  //
  G4Box* sWorld
    = new G4Box("World",                        //name
                 worldDimensions.x(),           //dimensions (half-lentghs)
		             worldDimensions.y(),
                 worldDimensions.z());

  G4LogicalVolume* worldVolume
    = new G4LogicalVolume(sWorld,               //shape
                          Vacuum,               //material
                          "World");             //name

  G4VPhysicalVolume* pWorld
    = new G4PVPlacement(0,                      //no rotation
                        G4ThreeVector(),        //at (0,0,0)
                        worldVolume,           //logical volume
                        "World",                //name
                        0,                      //mother  volume
                        false,                  //no boolean operation
                        0);                     //copy number

 
  worldVolume->SetMaterial(Vacuum);

  // "Collimator"
  if (strcmp( fCollType, "Rectangle" ) == 0 ) {
    //  Rectangular type: Make box and another box inside made of vacuum
    G4Box* coll
      = new G4Box("RColl",                        //its name
                 worldDimensions.x(),           //dimensions (half-lengths)
                 worldDimensions.y(),
		             worldDimensions.z());

    G4LogicalVolume* collVolume
      = new G4LogicalVolume(coll,                 //its shape
                          collMaterial,         //its material
                          "RColl");             //its name

    new G4PVPlacement(0,                         //no rotation
                    G4ThreeVector(),            //at (0,0,0)
                    collVolume,                 //its logical volume                           
                    "RColl",                    //its name
                    worldVolume,                //its mother  volume
                    false,                      //no boolean operation
                    0);                         //copy number

    G4Box* collInner
      = new G4Box("RCollInner",                   //its name
                 fCollAperX*m,                    //dimensions (half-lengths)
                 fCollAperY*m,
                 fCollLength*m);

    G4LogicalVolume* collVolumeInner
      = new G4LogicalVolume(collInner,            //its shape
                          Vacuum,               //its material
                          "RCollInner");        //its name

    new G4PVPlacement(0,                          //no rotation
                    G4ThreeVector(),            //at (0,0,0)
                    collVolumeInner,                //its logical volume                           
                    "RCollInner",                      //its name
                    collVolume,               //its mother  volume
                    false,                      //no boolean operation
                    0);                         //copy number  
  }
  else if (strcmp( fCollType, "Ellipse" ) ) {
    // Elliptical tube: surface equation: 1 = (x/Dx)^2 + (y/Dy)^2
    G4EllipticalTube* coll = new G4EllipticalTube( "EColl",                  // name
				 worldDimensions.x(),    // Dx
				 worldDimensions.y(),    // Dy
				 worldDimensions.z() );  // Dz
    G4LogicalVolume* collVolume
      = new G4LogicalVolume(coll,                 //its shape
                          collMaterial,         //its material
                          "EColl");             //its name
    new G4PVPlacement(0,                         //no rotation
                    G4ThreeVector(),            //at (0,0,0)
                    collVolume,                 //its logical volume                           
                    "EColl",                    //its name
                    worldVolume,                //its mother  volume
                    false,                      //no boolean operation
                    0);                         //copy number
    G4EllipticalTube* collInner = new G4EllipticalTube( "ECollInner",                  // name
				 fCollAperX*m,    // Dx
				 fCollAperY*m,    // Dy
				 fCollLength*m );  // Dz
    G4LogicalVolume* collVolumeInner
      = new G4LogicalVolume(collInner,            //its shape
                          Vacuum,               //its material
                          "ECollInner");        //its name

    new G4PVPlacement(0,                          //no rotation
                    G4ThreeVector(),            //at (0,0,0)
                    collVolumeInner,                //its logical volume                           
                    "ECollInner",                      //its name
                    collVolume,               //its mother  volume
                    false,                      //no boolean operation
                    0);                         //copy number  
  }
  else {
    G4cerr << "Unknown Geometry type requested" ;
    return NULL;
  }
  //always return the root volume
  //
  return pWorld;

}

