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
#include <string>

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
        fCollDY(lman->CollDY),
        fVacuumMaterial(lman->VacuumMaterial)
{
  // Define materials via NIST manager
  //
  nistManager = G4NistManager::Instance();
  
  // Store lucretiaManager class reference for field setup
  fLman=lman;
}

uint32_t geomConstruction::hashL(char *key, size_t len)
{
    uint32_t hash, i;
    for(hash = i = 0; i < len; ++i)
    {
        hash += key[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
}

G4Material* geomConstruction::ProcessMaterial(G4String materialName, int matID = 1 )
{
  // Generate GEANT4 material definition from given data
  static const double STP_Temperature = 273.15*kelvin;
  static const double STP_Pressure    = 1.*atmosphere;
  static std::vector<G4Material*> materialVector ;
  G4double density,a,z;
  G4double temperature, pressure, fractionmass;
  G4String name, symbol;
  static G4Material* vacuumMaterial=nistManager->ConstructNewGasMaterial("Vacuum","G4_H",273*kelvin,1.e-25*pascal);
  G4int ncomponents;
  G4Material* material ;
  //string eleChar ;
  //stringstream convert;
  //convert << fLman->fEle ;
  //eleChar = convert.str();
  if (strcmp(materialName,"Vacuum")==0) { // default perfect vacuum
    material = vacuumMaterial ;
  } // Process user generated material description(s)
  else if (strncmp(materialName,"User",4) == 0) {
    int imat = (int)(materialName[4]-'0') ;
    int ic ;
    density=fLman->UserMaterial[imat].density*g/cm3;
    if (fLman->UserMaterial[imat].pressure>0)
      pressure=fLman->UserMaterial[imat].pressure*pascal;
    else
      pressure=STP_Pressure;
    if (fLman->UserMaterial[imat].temperature>0)
      temperature=fLman->UserMaterial[imat].temperature*kelvin;
    else
      temperature=STP_Temperature;
    if ( strcmp(fLman->UserMaterial[imat].state,"Solid") == 0 )
      material = new G4Material(name=materialName, density, ncomponents=fLman->UserMaterial[imat].NumComponents, kStateSolid, temperature, pressure);
    else if ( strcmp(fLman->UserMaterial[imat].state,"Liquid") == 0 )
      material = new G4Material(name=materialName, density, ncomponents=fLman->UserMaterial[imat].NumComponents, kStateLiquid, temperature, pressure);
    else
      material = new G4Material(name=materialName, density, ncomponents=fLman->UserMaterial[imat].NumComponents, kStateGas, temperature, pressure);
    for (ic=0;ic<fLman->UserMaterial[imat].NumComponents;ic++)
      material->AddElement(new G4Element(name=fLman->UserMaterial[imat].element[ic].Name,
              symbol=fLman->UserMaterial[imat].element[ic].Symbol,
              z=fLman->UserMaterial[imat].element[ic].Z,
              fLman->UserMaterial[imat].element[ic].A*g/mole),
              fractionmass=fLman->UserMaterial[imat].element[ic].FractionMass) ;
  }
  else { // assume want something from the GEANT4 (NIST) material database
    if (fLman->MatP[matID]>0)
      pressure=fLman->MatP[matID]*pascal;
    else
      pressure=STP_Pressure;
    if (fLman->MatT[matID]>0)
      temperature=fLman->MatT[matID]*kelvin;
    else
      temperature=STP_Temperature;
    if (pressure!=STP_Pressure || temperature!=STP_Temperature)
      material = nistManager->ConstructNewGasMaterial(materialName,materialName,temperature,pressure);
    else
      material = nistManager->FindOrBuildMaterial(materialName);
  }
  return material ;
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
  fVacuumMaterial=lman->VacuumMaterial ;
  // Store lucretiaManager class reference for field setup
  fLman=lman;
}

geomConstruction::~geomConstruction()
{
}

G4VPhysicalVolume* geomConstruction::Construct()
{
  // Vacuum definition
  Vacuum = ProcessMaterial(fVacuumMaterial) ;
  // Collimator material
  collMaterial = ProcessMaterial(fCollMaterialName,1);
  collMaterial2 = ProcessMaterial(fCollMaterialName2,2);
  
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
    if (fCollAperX>0 && fCollAperY>0) {
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
    if (fCollAperX>0 && fCollAperY>0) {
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
  //if (fEmFieldSetup)
  //  delete fEmFieldSetup ;
  if (fLman->EnableEM==1) {
    if (fLman->Verbose>=1)
      printf("Generating EM field map...\n") ;
    FieldSetup* fEmFieldSetup = new FieldSetup(fLman);
  }
}


