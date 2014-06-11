#include "geomConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include <iostream>

geomConstruction::geomConstruction(lucretiaManager* lman, G4double dz)
: G4VUserDetectorConstruction(),
        fCollType(lman->GeomType),
        fCollMaterialName(lman->Material),
        fCollMaterialName2(lman->Material2),
        fCollAperX(lman->AperX),
        fCollAperY(lman->AperY),
        fCollThickness(lman->Thickness),
        fCollLength(dz),
        fCollLength2(lman->CollLen2),
        fCollAperX2(lman->AperX2),
        fCollAperY2(lman->AperY2),
        fCollAperX3(lman->AperX3),
        fCollAperY3(lman->AperY3),
        fCollDX(lman->CollDX),
        fCollDY(lman->CollDY)
{
  // Define materials via NIST manager
  //
  nistManager = G4NistManager::Instance();
  
  // Make Vacuum material for World volume
  G4double density,a,z;
  G4double temperature, pressure, fractionmass;
  G4String name, symbol;
  G4int ncomponents;
  density = 1.e-25*g/cm3;
  pressure    = 1.e-25*pascal;
  temperature = 2.73*kelvin;
  /*G4Material* Vacuum =
          new G4Material(name="Vacuum",      //Name as String
          1,		       //Atomic Number,  in this case we use 1 for hydrogen
          1.008*g/mole,  //Mass per Mole "Atomic Weight"  1.008*g/mole for Hydoren
          1.e-25*g/cm3,  //Density of Vaccuum  *Cant be Zero, Must be small insted
          kStateGas,     //kStateGas for Gas
          2.73*kelvin,   //Temperatuer for ga
          1.e-25*g/cm3); //Pressure for Vaccum */
  //a = 14.01*g/mole;
  //elN  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);
  //a = 16.00*g/mole;
  //elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);
  a = 1.008*g/mole;
  elH = new G4Element(name="Hydrogen"  ,symbol="H" , z= 1., a);
  Vacuum = new G4Material(name="Vacuum", density, ncomponents=1, kStateGas, temperature, pressure);
  //Vacuum->AddElement(elN, fractionmass=0.7);
  //Vacuum->AddElement(elO, fractionmass=0.3);
  Vacuum->AddElement(elH, fractionmass=1.);
  
  // Store lucretiaManager class reference for field setup
  fLman=lman;
}

void geomConstruction::SetGeomParameters(lucretiaManager* lman)
{
  fCollType=lman->GeomType;
  fCollMaterialName=lman->Material;
  fCollMaterialName2=lman->Material2;
  fCollAperX=lman->AperX;
  fCollAperY=lman->AperY;
  fCollThickness=lman->Thickness;
  fCollLength=lman->Lcut ;
  fCollLength2=lman->CollLen2 ;
  fCollAperX2=lman->AperX2 ;
  fCollAperY2=lman->AperY2 ;
  fCollAperX3=lman->AperX3 ;
  fCollAperY3=lman->AperY3 ;
  fCollDX=lman->CollDX ;
  fCollDY=lman->CollDY ;
}

geomConstruction::~geomConstruction()
{
}

G4VPhysicalVolume* geomConstruction::Construct()
{
  // Collimator material
  G4Material* collMaterial
          = nistManager->FindOrBuildMaterial(fCollMaterialName);
  G4Material* collMaterial2
          = nistManager->FindOrBuildMaterial(fCollMaterialName2);
  
  // Geometry parameters
  //
  G4ThreeVector worldDimensions = G4ThreeVector((fCollAperX+fCollThickness)*m,(fCollAperY+fCollThickness)*m,fCollLength*m);
  
  // Clean old geometry, if any //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  // World
  sWorld
          = new G4Box("World",                        //name
          worldDimensions.x(),           //dimensions (half-lentghs)
          worldDimensions.y(),
          worldDimensions.z());
  
  worldVolume
          = new G4LogicalVolume(sWorld,               //shape
          Vacuum,               //material
          "World");             //name
  
  pWorld
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
    collBox
            = new G4Box("RColl",                        //its name
            worldDimensions.x(),           //dimensions (half-lengths)
            worldDimensions.y(),
            worldDimensions.z());
    
    collVolume
            = new G4LogicalVolume(collBox,                 //its shape
            collMaterial,         //its material
            "RColl");             //its name
    
    collPlacement
            = new G4PVPlacement(0,                         //no rotation
            G4ThreeVector(),            //at (0,0,0)
            collVolume,                 //its logical volume
            "RColl",                    //its name
            worldVolume,                //its mother  volume
            false,                      //no boolean operation
            0);                         //copy number
    
    collInnerBox
            = new G4Box("RCollInner",                   //its name
            fCollAperX*m,                    //dimensions (half-lengths)
            fCollAperY*m,
            fCollLength*m);
    
    collVolumeInner
            = new G4LogicalVolume(collInnerBox,            //its shape
            Vacuum,               //its material
            "RCollInner");        //its name
    collInnerPlacement
            = new G4PVPlacement(0,                          //no rotation
            G4ThreeVector(),            //at (0,0,0)
            collVolumeInner,                //its logical volume
            "RCollInner",                      //its name
            collVolume,               //its mother  volume
            false,                      //no boolean operation
            0);                         //copy number
  }
  else if (strcmp( fCollType, "Ellipse" ) == 0) {
    // Elliptical tube: surface equation: 1 = (x/Dx)^2 + (y/Dy)^2
    collTube = new G4EllipticalTube( "EColl",                  // name
            worldDimensions.x(),    // Dx
            worldDimensions.y(),    // Dy
            worldDimensions.z() );  // Dz
    collVolume
            = new G4LogicalVolume(collTube,                 //its shape
            collMaterial,         //its material
            "EColl");             //its name
    collPlacement
            = new G4PVPlacement(0,                         //no rotation
            G4ThreeVector(),            //at (0,0,0)
            collVolume,                 //its logical volume
            "EColl",                    //its name
            worldVolume,                //its mother  volume
            false,                      //no boolean operation
            0);                         //copy number
    collInnerTube = new G4EllipticalTube( "ECollInner",                  // name
            fCollAperX*m,    // Dx
            fCollAperY*m,    // Dy
            fCollLength*m );  // Dz
    collVolumeInner
            = new G4LogicalVolume(collInnerTube,            //its shape
            Vacuum,               //its material
            "ECollInner");        //its name
    collInnerPlacement
            = new G4PVPlacement(0,                          //no rotation
            G4ThreeVector(),            //at (0,0,0)
            collVolumeInner,                //its logical volume
            "ECollInner",                      //its name
            collVolume,               //its mother  volume
            false,                      //no boolean operation
            0);                         //copy number
  }
  else if (strcmp( fCollType, "Tapered" ) == 0 ) {
    //  Tapered Collimator with optional box insert
    collTapBox
            = new G4Box("TapColl",                        //its name
            worldDimensions.x(),           //dimensions (half-lengths)
            worldDimensions.y(),
            worldDimensions.z());
    
    collVolume
            = new G4LogicalVolume(collTapBox,                 //its shape
            collMaterial,         //its material
            "TapColl");             //its name
    
    collPlacement
            = new G4PVPlacement(0,                         //no rotation
            G4ThreeVector(),            //at (0,0,0)
            collVolume,                 //its logical volume
            "TapColl",                    //its name
            worldVolume,                //its mother  volume
            false,                      //no boolean operation
            0);                         //copy number
    double tapLength = (fCollLength-fCollLength*fCollLength2) ;
    collTapBox1
            = new G4Trd("TapColl1",
               (fCollAperX+fCollAperX2*fCollDX)*m,
               fCollAperX*m,
               (fCollAperY+fCollAperY2*fCollDY)*m,
               fCollAperY*m,
               (tapLength/2)*m );
    if (fCollLength2>0) {
      collTapBox2
              = new G4Box("TapColl2",
                 (fCollAperX+fCollAperX3*fCollDX)*m,
                 (fCollAperY+fCollAperY3*fCollDY)*m,
                 (fCollLength*fCollLength2)*m);
      collTapBox3
              = new G4Box("TapColl3",
                 fCollAperX*m,
                 fCollAperY*m,
                 (fCollLength*fCollLength2)*m);
    }
    collTapBox4
            = new G4Trd("TapColl4",
               fCollAperX*m,
               (fCollAperX+fCollAperX2*fCollDX)*m,
               fCollAperY*m,
               (fCollAperY+fCollAperY2*fCollDY)*m,
               (tapLength/2)*m);
    collInnerPlacement
            = new G4PVPlacement(0,    //no rotation
            G4ThreeVector(0,0,(-fCollLength+tapLength/2)*m),
            new G4LogicalVolume(collTapBox1,Vacuum,"TapCollInner1"),
            "TapCollInner1",                      //its name
            collVolume,               //its mother  volume
            false,                      //no boolean operation
            0);      //copy number
    collInnerPlacement2
            = new G4PVPlacement(0,    //no rotation
            G4ThreeVector(0,0,(fCollLength*fCollLength2+tapLength/2)*m),
            new G4LogicalVolume(collTapBox4,Vacuum,"TapCollInner2"),
            "TapCollInner2",                      //its name
            collVolume,               //its mother  volume
            false,                      //no boolean operation
            0);                         //copy number
    if (fCollLength2>0) {
      collVolumeInner
              = new G4LogicalVolume(collTapBox2,            //its shape
              collMaterial2,               //its material
              "TapCollInner3");        //its name
      collInnerPlacement3
            = new G4PVPlacement(0,                          //no rotation
            G4ThreeVector(),            //at (0,0,0)
            collVolumeInner,                //its logical volume
            "TapCollInner3",                      //its name
            collVolume,               //its mother  volume
            false,                      //no boolean operation
            0);                         //copy number
      collInnerPlacement4
            = new G4PVPlacement(0,                          //no rotation
            G4ThreeVector(),            //at (0,0,0)
            new G4LogicalVolume(collTapBox3,Vacuum,"TapCollInner4"),
            "TapCollInner4",                      //its name
            collVolumeInner,               //its mother  volume
            false,                      //no boolean operation
            0);                         //copy number
    }
    
  }
  else {
    G4cerr << "Unknown Geometry type requested" ;
    return NULL;
  }
  
  //always return the root volume
  //
  return pWorld;
  
}
void geomConstruction::ConstructSDandField()
{
  // Construct the field creator - this will register the field it creates
  if (fEmFieldSetup)
    delete fEmFieldSetup ;
  if (fLman->EnableEM==1) {
    if (fLman->Verbose>=1)
      printf("Generating EM field map...\n") ;
    FieldSetup* fEmFieldSetup = new FieldSetup(fLman);
  }
}


