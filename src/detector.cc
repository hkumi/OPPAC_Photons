#include "detector.hh"


MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{

}

MySensitiveDetector::~MySensitiveDetector()
{}

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{

    G4Track *track = aStep->GetTrack();

    //track->SetTrackStatus(fStopAndKill);
    G4ParticleDefinition*  particle =aStep->GetTrack()->GetDefinition();
    if (particle->GetParticleName() == "opticalphoton"){   
       G4StepPoint *preStepPoint = aStep->GetPreStepPoint();//used when the neutron enters the detector
       G4StepPoint *postStepPoint = aStep->GetPostStepPoint();//used when the neutron leaves the detector
       const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();

       G4ThreeVector posPhotons = preStepPoint->GetPosition();//accessing the position       
       G4String particleName = particle->GetParticleName();
       G4cout << "This is the particle in the sensor: " << particleName << G4endl;       


       G4int copyNo = touchable->GetCopyNumber();
       G4cout << "Copy number: " << copyNo << G4endl;

       G4VPhysicalVolume *physVol = touchable->GetVolume();
       G4ThreeVector posDetector = physVol->GetTranslation();
       G4cout<< "Neutron Position:"<< posPhotons << G4endl;

       G4AnalysisManager *man = G4AnalysisManager::Instance();

       man->FillNtupleDColumn(0, 0, posPhotons[0]/mm);
       man->AddNtupleRow(0);

       man->FillNtupleDColumn(1, 0, posPhotons[1]/mm);
  
       man->AddNtupleRow(1);

       man->FillNtupleDColumn(2, 0, copyNo);
       man->AddNtupleRow(2);
 

       man->FillH2(0, posPhotons[0]/mm, posPhotons[1]/mm);


   }
/*
       man->FillNtupleDColumn(2, 0, wave);
       man->AddNtupleRow(2);
 
       man->FillNtupleDColumn(3, 0, tg);
       man->AddNtupleRow(3);


       man->FillH2(0, X/mm, Y/mm);

*/

    //}
/*
    G4StepPoint *preStepPoint = aStep->GetPreStepPoint();//used when the neutron enters the detector
    G4StepPoint *postStepPoint = aStep->GetPostStepPoint();//used when the neutron leaves the detector

    G4ThreeVector posNeutron = preStepPoint->GetPosition();//accessing the position
    G4ThreeVector momNeutron = preStepPoint->GetMomentum();//accessing the momentum. 
    //G4String energyString = G4BestUnit(energy, "Energy"); // Assuming you want the energy in the best unit for energy
    G4double wlen = (1.239841939*eV/momNeutron.mag())*1E+03; //the deBroglie wavelength of a neutron in nm

    const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();

    G4int copyNo = touchable->GetCopyNumber();

    //G4cout << "Copy number: " << copyNo << G4endl;

    G4VPhysicalVolume *physVol = touchable->GetVolume();
    G4ThreeVector posDetector = physVol->GetTranslation();

    G4cout<< "Neutron Position:"<< posDetector << G4endl;

    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

    G4AnalysisManager *man = G4AnalysisManager::Instance();
    man->FillNtupleIColumn(0, 0, evt);
    man->FillNtupleDColumn(0, 1, posNeutron[0]);
    man->FillNtupleDColumn(0, 2, posNeutron[1]);
    man->FillNtupleDColumn(0, 3, posNeutron[2]);
    man->FillNtupleDColumn(0, 4, wlen);
    man->AddNtupleRow(0);

    man->FillNtupleIColumn(1, 0, evt);
    man->FillNtupleDColumn(1, 1, posDetector[0]);
    man->FillNtupleDColumn(1, 2, posDetector[1]);
    man->FillNtupleDColumn(1, 3, posDetector[2]);
    man->AddNtupleRow(1);
*/


 
    //G4ThreeVector momNeutron = preStepPoint->GetMomentum();
   /*
    G4double time = preStepPoint->GetGlobalTime();

    G4double wlen = (1.239841939*eV/momPhoton.mag())*1E+03;


    #ifndef G4MULTITHREADED
        //G4cout << "Detector position: " << posDetector << G4endl;
        //G4cout << "Photon wavelength: " << wlen << G4endl;
    #endif*/

 return true;
}
