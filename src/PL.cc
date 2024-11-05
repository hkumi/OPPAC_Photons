#include "globals.hh"
#include "PL.hh"
#include "G4SystemOfUnits.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PL::PL():  G4VUserPhysicsList()
{
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PL::~PL()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PL::ConstructParticle()
{

 // gamma
  G4Gamma::Gamma();

  // leptons
  G4Electron::Electron();
  G4Positron::Positron();

  // ions
  G4Deuteron::Deuteron();
  G4Triton::Triton();
  G4He3::He3();
  G4Alpha::Alpha();
  G4GenericIon::GenericIonDefinition();

  // optical photons
  G4OpticalPhoton::OpticalPhotonDefinition();
}


void PL::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructOptical();
  ConstructGeneral();
  AddStepMax();
  SetCuts();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4ionIonisation.hh"

#include "G4Scintillation.hh"
#include "G4Cerenkov.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpWLS.hh"

#include "G4IonParametrisedLossModel.hh"
//#include "G4EmProcessOptions.hh"
#include "G4EmParameters.hh"
#include "G4LossTableManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void PL::ConstructEM()
{
  // Use G4EmParameters instead of G4EmProcessOptions
  G4EmParameters* emParams = G4EmParameters::Instance();
/*  emParams->SetVerbose(0);
  emParams->SetMaxEnergy(100 * GeV);
  emParams->SetNumberOfBinsPerDecade(200);
  emParams->SetLambdaBinning(200);
  emParams->SetFluo(true);   // Enable fluorescence
  emParams->SetPIXE(true);   // Enable PIXE (Particle-Induced X-ray Emission)
  emParams->SetAuger(true);  // Enable Auger effect
  emParams->SetPixe(true); // Enable PIXE, if SetPixe exists in your version

  auto theParticleIterator = GetParticleIterator();
  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    pmanager->SetVerboseLevel(0);
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {
      // gamma
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4ComptonScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversion);

    } else if (particleName == "e-") {
      // electron
      pmanager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4eIonisation, -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung, -1, 3, 3);

    } else if (particleName == "e+") {
      // positron
      pmanager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4eIonisation, -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung, -1, 3, 3);
      pmanager->AddProcess(new G4eplusAnnihilation, 0, -1, 4);

    } else if (particleName == "alpha" || particleName == "deuteron" ||
               particleName == "triton" || particleName == "He3") {
      // alpha, deuteron, triton, and He3
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      pmanager->AddProcess(msc, -1, 1, 1);

      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetEmModel(new G4IonParametrisedLossModel);
      ionIoni->SetStepFunction(0.1, 20 * um);
      pmanager->AddProcess(ionIoni, -1, 2, 2);

    } else if (particleName == "GenericIon") {
      // genericIon
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);

      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetEmModel(new G4IonParametrisedLossModel());
      ionIoni->SetStepFunction(0.1, 20 * um);
      pmanager->AddProcess(ionIoni, -1, 2, 2);
    }
  }*/
}

// Optical Processes ////////////////////////////////////////////////////////
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
//#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"

void PL::ConstructOptical() 
{

  auto theParticleIterator = GetParticleIterator();
/*  // default scintillation process
  G4Scintillation* theScintProcessDef = new G4Scintillation("Scintillation");
  // theScintProcessDef->DumpPhysicsTable();
  theScintProcessDef->SetTrackSecondariesFirst(true);
  theScintProcessDef->SetScintillationYieldFactor(1.0); //
  theScintProcessDef->SetScintillationExcitationRatio(0.0); //
  theScintProcessDef->SetVerboseLevel(0);

  // scintillation process for alpha:
  G4Scintillation* theScintProcessAlpha = new G4Scintillation("Scintillation");
  // theScintProcessNuc->DumpPhysicsTable();
  theScintProcessAlpha->SetTrackSecondariesFirst(true);
  theScintProcessAlpha->SetScintillationYieldFactor(1.1);
  theScintProcessAlpha->SetScintillationExcitationRatio(1.0);
  theScintProcessAlpha->SetVerboseLevel(0);

  // scintillation process for heavy nuclei
  G4Scintillation* theScintProcessNuc = new G4Scintillation("Scintillation");
  // theScintProcessNuc->DumpPhysicsTable();
  theScintProcessNuc->SetTrackSecondariesFirst(true);
  theScintProcessNuc->SetScintillationYieldFactor(0.2);
  theScintProcessNuc->SetScintillationExcitationRatio(1.0);
  theScintProcessNuc->SetVerboseLevel(0);

  // optical processes
  G4OpAbsorption* theAbsorptionProcess = new G4OpAbsorption();
  //  G4OpRayleigh* theRayleighScatteringProcess = new G4OpRayleigh();
  G4OpBoundaryProcess* theBoundaryProcess = new G4OpBoundaryProcess();
  //  theAbsorptionProcess->DumpPhysicsTable();
  //  theRayleighScatteringProcess->DumpPhysicsTable();
  theAbsorptionProcess->SetVerboseLevel(0);
  // theRayleighScatteringProcess->SetVerboseLevel(OpVerbLevel);
  theBoundaryProcess->SetVerboseLevel(0);

  theParticleIterator->reset();
  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
      if (theScintProcessDef->IsApplicable(*particle)) {
	//      if(particle->GetPDGMass() > 5.0*GeV) 
	if(particle->GetParticleName() == "GenericIon") {
	  pmanager->AddProcess(theScintProcessNuc); // AtRestDiscrete
	  pmanager->SetProcessOrderingToLast(theScintProcessNuc,idxAtRest);
	  pmanager->SetProcessOrderingToLast(theScintProcessNuc,idxPostStep);
	}	  
	else if(particle->GetParticleName() == "alpha") {
	  pmanager->AddProcess(theScintProcessAlpha);
	  pmanager->SetProcessOrderingToLast(theScintProcessAlpha,idxAtRest);
	  pmanager->SetProcessOrderingToLast(theScintProcessAlpha,idxPostStep);
	}
	else {
	  pmanager->AddProcess(theScintProcessDef);
	  pmanager->SetProcessOrderingToLast(theScintProcessDef,idxAtRest);
	  pmanager->SetProcessOrderingToLast(theScintProcessDef,idxPostStep);
	}	  
      }
      
      if (particleName == "opticalphoton") {
	pmanager->AddDiscreteProcess(theAbsorptionProcess);
	//	pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
	pmanager->AddDiscreteProcess(theBoundaryProcess);
      }
    }*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Decay.hh"

void PL::ConstructGeneral()
{
  auto theParticleIterator = GetParticleIterator();
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    pmanager->SetVerboseLevel(0);
    if (theDecayProcess->IsApplicable(*particle)) {
      pmanager ->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"

void PL::AddStepMax()
{
  auto theParticleIterator = GetParticleIterator();
  // Step limitation seen as a process
  G4StepLimiter* stepLimiter = new G4StepLimiter();
  theParticleIterator->reset();
  while ((*theParticleIterator)())
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    pmanager->SetVerboseLevel(0);
    if (particle->GetPDGCharge() != 0.0)
    {
      pmanager ->AddDiscreteProcess(stepLimiter);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PL::SetCuts()
{
  //G4VUserPL::SetCutsWithDefault method sets
  //the default cut value for all particle types
  //
  SetCutsWithDefault();
  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
