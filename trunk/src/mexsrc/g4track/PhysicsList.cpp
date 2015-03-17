#include "PhysicsList.hh"
//#include "PhysicsListMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4DecayPhysics.hh"
#include "G4Decay.hh"
#include "G4DecayTable.hh"
#include "G4DecayWithSpin.hh"
#include "G4MuonDecayChannelWithSpin.hh"
#include "G4MuonRadiativeDecayChannelWithSpin.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4GammaConversionToMuons.hh"
#include "G4AnnihiToMuPair.hh"
#include "G4StoppingPhysics.hh"
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4ElectronNuclearProcess.hh"
#include "G4MuonNuclearProcess.hh"
#include "G4PhotoNuclearProcess.hh"
#include "G4PositronNuclearProcess.hh"
#include "G4eeToHadrons.hh"

#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4SynchrotronRadiation.hh"
#include "G4SynchrotronRadiationInMat.hh"

#include "G4EmProcessOptions.hh"

#include "G4StepLimiter.hh"

#include "G4SystemOfUnits.hh"

PhysicsList::PhysicsList(): G4VModularPhysicsList()
{
  //fMess = new PhysicsListMessenger(this);
  defaultCutValue = 1.*km;
  fSRType = true; 
  fParticleList = new G4DecayPhysics("decays") ;
  RegisterPhysics( new G4HadronPhysicsFTFP_BERT() ) ;
}

PhysicsList::~PhysicsList()
{ 
  //delete fMess;
  delete fParticleList;
}

void PhysicsList::ConstructParticle()
{
  fParticleList->ConstructParticle();
  /*
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();

  
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();

  
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  
  G4Gamma::GammaDefinition();
  G4OpticalPhoton::OpticalPhotonDefinition();
  */

}

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  theParticleIterator->reset();
  G4EmProcessOptions opt;
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4ComptonScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversion);
      pmanager->AddDiscreteProcess(new G4RayleighScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversionToMuons);
      pmanager->AddDiscreteProcess(new G4PhotoNuclearProcess) ;
    } else if (particleName == "e-") {
      pmanager->AddDiscreteProcess(new G4ElectronNuclearProcess) ;
      pmanager->AddProcess(new G4eMultipleScattering,       -1, 1, 1);
      pmanager->AddProcess(new G4eIonisation,               -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung,           -1, 3, 3);
      if (fSRType) {
	pmanager->AddProcess(new G4SynchrotronRadiation,      -1,-1, 4);
      } else {
	pmanager->AddProcess(new G4SynchrotronRadiationInMat, -1,-1, 4); 
      }
      pmanager->AddProcess(new G4StepLimiter,               -1,-1, 5);
    } else if (particleName == "e+") {
      pmanager->AddProcess(new G4eMultipleScattering,       -1, 1, 1);
      pmanager->AddProcess(new G4eIonisation,               -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung,           -1, 3, 3);
      pmanager->AddProcess(new G4eplusAnnihilation,          0,-1, 4);
      if (fSRType) {
	pmanager->AddProcess(new G4SynchrotronRadiation,      -1,-1, 5);
      } else {
	pmanager->AddProcess(new G4SynchrotronRadiationInMat, -1,-1, 5);
      }
      pmanager->AddProcess(new G4StepLimiter,               -1,-1, 6);
      pmanager->AddDiscreteProcess(new G4AnnihiToMuPair);
      pmanager->AddDiscreteProcess(new G4PositronNuclearProcess);
      pmanager->AddDiscreteProcess(new G4eeToHadrons);
    } else if( particleName == "mu+" || particleName == "mu-"    ) {
      pmanager->AddDiscreteProcess(new G4MuonNuclearProcess) ;
      pmanager->AddProcess(new G4MuMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4MuIonisation,         -1, 2, 2);
      pmanager->AddProcess(new G4MuBremsstrahlung,     -1, 3, 3);
      pmanager->AddProcess(new G4MuPairProduction,     -1, 4, 4);
    }
    opt.SetVerbose(0);
  }
  ConstructGeneral();
}

#include "G4Decay.hh"

void PhysicsList::ConstructGeneral()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while ((*theParticleIterator)()){
      G4ParticleDefinition* particle = theParticleIterator->value();      G4ProcessManager* pmanager = particle->GetProcessManager();
      if (theDecayProcess->IsApplicable(*particle)) {
        pmanager ->AddProcess(theDecayProcess);
        // set ordering for PostStepDoIt and AtRestDoIt
        pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
        pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
      }
  }
}


void PhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");

  if (verboseLevel>0) DumpCutValuesTable();
}




