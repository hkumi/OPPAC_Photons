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
       //G4cout << "This is the particle in the sensor: " << particleName << G4endl;       


       G4int copyNo = touchable->GetCopyNumber();

       G4String volumeName = touchable->GetVolume()->GetName(); // name of the volume
       //G4cout << "Photon detected in " << volumeName << "with Copy Number" << copyNo<< G4endl;
       G4AnalysisManager *man = G4AnalysisManager::Instance();


       // taking photons in the different sensor arrays based on the volume name
       if (volumeName == "sensor_Vol1") {
            // for sensor_Vol1
          G4cout << "Photon detected in sensor_Vol1 with Copy Number: " << copyNo << G4endl;
          man->FillNtupleDColumn(0, 0, posPhotons[0]/mm);
          man->FillNtupleDColumn(0, 1, posPhotons[1]/mm);
          man->FillNtupleDColumn(0, 2, posPhotons[2]/mm);
          man->FillNtupleIColumn(0, 3, copyNo);
          man->AddNtupleRow(0);

          } else if (volumeName == "sensor_Vol2") {
            // for sensor_Vol2
            G4cout << "Photon detected in sensor_Vol2 with Copy Number: " << copyNo << G4endl;
            man->FillNtupleDColumn(1, 0, posPhotons[0]/mm);
            man->FillNtupleDColumn(1, 1, posPhotons[1]/mm);
            man->FillNtupleDColumn(1, 2, posPhotons[2]/mm);
            man->FillNtupleIColumn(1, 3, copyNo);
            man->AddNtupleRow(1);

        } else if (volumeName == "sensor_Vol3") {
            //  sensor_Vol3
            G4cout << "Photon detected in sensor_Vol3 with Copy Number: " << copyNo << G4endl;
            man->FillNtupleDColumn(2, 0, posPhotons[0]/mm);
            man->FillNtupleDColumn(2, 1, posPhotons[1]/mm);
            man->FillNtupleDColumn(2, 2, posPhotons[2]/mm);
            man->FillNtupleIColumn(2, 3, copyNo);
            man->AddNtupleRow(2);

        } else if (volumeName == "sensor_Vol4") {
            // sensor_Vol3
            G4cout << "Photon detected in sensor_Vol4 with Copy Number: " << copyNo << G4endl;
            man->FillNtupleDColumn(3, 0, posPhotons[0]/mm);
            man->FillNtupleDColumn(3, 1, posPhotons[1]/mm);
            man->FillNtupleDColumn(3, 2, posPhotons[2]/mm);
            man->FillNtupleIColumn(3, 3, copyNo);
            man->AddNtupleRow(3);
       
          
       }     


    }


 return true;
}
