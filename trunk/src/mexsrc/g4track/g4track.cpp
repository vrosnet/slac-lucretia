//extern "C" {
#include "../LucretiaCommon.h"
//}
#include "g4track.h"
#ifndef LUCRETIA_MANAGER
#include "lucretiaManager.hh"
#endif
//#include "G4MTRunManager.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "geomConstruction.hh"
#include "actionInitialization.hh"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include <iostream>

int g4track(int* blele, struct Bunch* ThisBunch, double* L)
{
  static int runno=1;
  // lucretiaManager manages interaction with Lucretia bunch structure and interfaces with Matlab data structures
  lucretiaManager* lman = new lucretiaManager(blele, ThisBunch) ;
  // return with error if no associated EXT Process on this BEAMLINE or there is a problem with one of the EXT Process object's properties
  if (lman->Status!=0)
    return lman->Status ;

  // construct the default run manager (the first time this function is called only)
  //static G4MTRunManager* runManager;
  static G4RunManager* runManager ;
  if (runManager == NULL) {
    runManager = new G4RunManager;
  
    // Collimator geometry definition (apertures and lengths are half-lengths in meters)
    G4double length = *L ;
    // Cuts to determine which particles get returned to Lucretia tracking after geant tracking through the geometry
    G4double ecut = lman->Ecut; // GeV
    G4double zcut = length ; // m
    // collType ("Rectangle" or "Ellipse"), materialName, aper_x, aper_y, thickness, length
    geomConstruction* thisGeomConstruction = new geomConstruction(lman->GeomType, lman->Material, lman->AperX, lman->AperY,
                  lman->Thickness, length);
    if (thisGeomConstruction == NULL)
      return -2 ;
    runManager->SetUserInitialization(thisGeomConstruction);

    // Setup physics processes we wish to use
    G4VModularPhysicsList* physicsList = new FTFP_BERT; // Default ALL physics process list
    //physicsList->RegisterPhysics(new G4StepLimiterPhysics());
    runManager->SetUserInitialization(physicsList);

    // Initialize action routines (to set primary particles and routine to store tracking results)
    runManager->SetUserInitialization(new actionInitialization(length,ecut,zcut,lman)); //

    // initialize G4 kernel
    runManager->Initialize();

    // get the pointer to the UI manager and set verbosities
    G4UImanager* UI = G4UImanager::GetUIpointer();
 
    if (lman->Verbose>0) {
      UI->ApplyCommand("/run/verbose 1");
      UI->ApplyCommand("/event/verbose 1");
      UI->ApplyCommand("/tracking/verbose 1");
    }
    else {
      UI->ApplyCommand("/run/verbose 0");
      UI->ApplyCommand("/event/verbose 0");
      UI->ApplyCommand("/tracking/verbose 0");
    }
  }
  // start a run
  runManager->BeamOn(runno);
  runno++;

  // clean up and return
  //delete runManager;
  delete lman;
  return 0 ;
}

