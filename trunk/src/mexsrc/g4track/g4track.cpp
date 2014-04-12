#include "../LucretiaCommon.h"
#include "g4track.h"
#ifndef LUCRETIA_MANAGER
#include "lucretiaManager.hh"
#endif
//#include "G4MTRunManager.hh" // -> this for the multi-threaded version
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "geomConstruction.hh"
#include "actionInitialization.hh"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include <string>

using namespace std ;

// Interface from Lucretia to GEANT4 tracking routines
// - returns number of Lucretia macro-particles un-stopped
int g4track(int* blele, int* bunchno, struct Beam* TheBeam, double* L)
{                     
  struct Bunch* ThisBunch = TheBeam->bunches[*bunchno] ;
  static int runno=1; // run number starts at 1

  // construct the run manager, geometry, physics, UI and action objects (the first time this function is called only)
  //static G4MTRunManager* runManager; // -> this for the multi-threaded version
  static G4RunManager* runManager ;
  static G4UImanager* UI ;
  static geomConstruction* thisGeomConstruction ;
  static G4VModularPhysicsList* physicsList ;
  static actionInitialization* thisAction ;
  static lucretiaManager* lman ;
  static int firstCall=1;
  // Collimator geometry definition (apertures and lengths are half-lengths in meters)
  G4double length = *L ;
  if (runManager == NULL) {
    // lucretiaManager manages interaction with Lucretia bunch structure and interfaces with Matlab data structures
    lman = new lucretiaManager(blele, bunchno, ThisBunch, *L) ;
    // return with error if no associated EXT Process on this BEAMLINE or there is a problem with one of the EXT Process object's properties
    if (lman->Status!=0)
      return lman->Status ;
    runManager = new G4RunManager;
    // get the pointer to the UI manager
    UI = G4UImanager::GetUIpointer();
    // Setup physics processes we wish to use (for now this is just all of them)
    physicsList = new FTFP_BERT; // Default ALL physics process list
    runManager->SetUserInitialization(physicsList);
    thisGeomConstruction = new geomConstruction(lman->GeomType, lman->Material, lman->AperX, lman->AperY,
                lman->Thickness, length);
    if (thisGeomConstruction == NULL)
      return -2 ;
    runManager->SetUserInitialization(thisGeomConstruction);
    // Initialize action routines (to set primary particles and routine to store tracking results)
    thisAction = new actionInitialization(lman) ;
    runManager->SetUserInitialization(thisAction); //
  }
  else { // If this function already called once, just re-initialise the Lucretia pointers to the new element and geometry
    lman->Initialize(blele, bunchno, ThisBunch, *L) ;
    if (lman->Status!=0)
      return lman->Status ;
    runManager->ReinitializeGeometry(); // Force new run to set new geometry
    thisGeomConstruction->SetGeomParameters(lman) ; // Create the new geometry
  }
  
  // initialize G4 kernel
  if (firstCall==0) {
    //UI->ApplyCommand("/run/physicsModified");
  }
  else {
    firstCall = 0;
    runManager->Initialize();
  }
  // set verbosities
  if (lman->Verbose==1) {
    UI->ApplyCommand("/run/verbose 1");
    UI->ApplyCommand("/event/verbose 0");
    UI->ApplyCommand("/tracking/verbose 0");
    UI->ApplyCommand("/control/verbose 2");
    UI->ApplyCommand("/run/dumpRegion") ;
  }
  else if (lman->Verbose==2) {
    UI->ApplyCommand("/run/verbose 2");
    UI->ApplyCommand("/event/verbose 1");
    UI->ApplyCommand("/tracking/verbose 1");
    UI->ApplyCommand("/control/verbose 2");
    UI->ApplyCommand("/run/dumpRegion") ;
  }
  else {
    UI->ApplyCommand("/run/verbose 0");
    UI->ApplyCommand("/event/verbose 0");
    UI->ApplyCommand("/tracking/verbose 0");
    UI->ApplyCommand("/control/verbose 0");
  }
  // start a run
  runManager->BeamOn(runno);
  runno++;

  // return
  // - return variables to Matlab workspace
  lman->SetLucretiaData() ;
  return lman->fNumRaysResumed ;
}

