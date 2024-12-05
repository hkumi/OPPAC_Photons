#include "detector.hh"
#include <numeric>
#include <map>
#include <vector>
#include <fstream>  // For writing to a file
#include <sstream>  // For string stream if you need formatting
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include "TF1.h"
#include "G4SDManager.hh"
#include "G4String.hh"
#include "G4HCofThisEvent.hh"
#include "SensorHit.hh"

MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{
 G4String HCname = "SensorCollection";
 collectionName.insert(HCname);
}

MySensitiveDetector::~MySensitiveDetector()
{}


void MySensitiveDetector::Initialize(G4HCofThisEvent* HCE) {
    // Create the hits collection

    SensorCollection = new SensorHitsCollection(SensitiveDetectorName, collectionName[0]);
    G4cout << "SensitiveDetectorName: " << SensitiveDetectorName << G4endl;
G4cout << "Collection Name: " << collectionName[0] << G4endl;

    // Register the collection with the event
    static G4int HCID = -1; // Static to ensure it's set only once
    if (HCID < 0) {
        HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }
    if (HCE) {
        HCE->AddHitsCollection(HCID, SensorCollection);
    } else {
        G4cerr << "Error: HCE is null in Initialize!" << G4endl;
    }
}



// Define the function to save positions and copy numbers to a CSV file
void MySensitiveDetector::SaveToCSV(G4ThreeVector posPhotons, G4int copyNo) {
    // Open the CSV file in append mode (to keep adding rows)
    std::ofstream outputFile("photon_positions_4.csv", std::ios::app);  

    // Check if the file is open
    if (outputFile.is_open()) {
        // Write the data to the CSV file
        // Each entry is written in a new line with positions (x, y, z) and copy number
        outputFile << posPhotons.x()/mm << "," 
                   << posPhotons.y()/mm << "," 
                   << posPhotons.z()/mm << "," 
                   << copyNo << "\n";
    } else {
        G4cerr << "Failed to open file!" << G4endl;
    }

    // Close the file
    outputFile.close();
}

void MySensitiveDetector::RecordSensorData(const std::string& volumeName, int ntupleIndex, double posX, double posY, int event, int copyNo, G4AnalysisManager* man) {
    // Fill ntuple with data
    man->FillNtupleDColumn(ntupleIndex, 0, posX / mm);   // X position
    man->FillNtupleDColumn(ntupleIndex, 1, posY / mm);   // Y position
    man->FillNtupleDColumn(ntupleIndex, 2, event);       // Event number
    man->FillNtupleIColumn(ntupleIndex, 3, copyNo);      // Sensor copy number
    man->AddNtupleRow(ntupleIndex);
}

// Calculate the weighted mean for the x-coordinate
double  MySensitiveDetector::calculateWeightedMeanX(double Px1, double Nx1, double sigmax1) {
    if (sigmax1 == 0 ) {
        std::cout << "Error: Sigma values cannot be zero!" << G4endl;
        return -1; // Return a sentinel value to indicate error
    }
    return (Px1 * Nx1 / sigmax1 ) / (Nx1 / sigmax1 );
}



void MySensitiveDetector::FitHistogram(const std::vector<int>& copyNumbers) {
    if (copyNumbers.empty()) return;

    // Find min and max copy numbers
    double min = *std::min_element(copyNumbers.begin(), copyNumbers.end());
    double max = *std::max_element(copyNumbers.begin(), copyNumbers.end());
    int bins = 50;

    // Create histogram
    
     static int histCounter = 0; // Static counter for unique histogram names
    TH1D* hist = new TH1D(Form("copy_%d", histCounter++), "Copy Number Distribution", bins, min, max);
    for (int copyNo : copyNumbers) {
        hist->Fill(copyNo);
    }

    // Fit the histogram with a Gaussian
    TF1* gaussFit1 = new TF1("gaussFit", "gaus", min, max);
    gaussFit1->SetLineColor(kRed);
    hist->Fit(gaussFit1, "R");

    // Extract Gaussian fit parameters
    double mean = gaussFit1->GetParameter(1);      // Mean (μ)
    double sigma = gaussFit1->GetParameter(2);     // Standard deviation (σ)
    double fwhm = 2.355 * sigma;                   // Full Width at Half Maximum (FWHM)
    double amplitude = gaussFit1->GetParameter(0); // Amplitude (height of the peak)

    // Log results (or save to ntuple, etc.)
    G4cout << "Gaussian Fit Results: "
           << "Mean: " << mean << ", Sigma: " << sigma
           << ", FWHM: " << fwhm << ", Amplitude: " << amplitude << G4endl;

    // Calculate x_position
    double x_pos = calculateWeightedMeanX( mean, amplitude, sigma);
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    man->FillNtupleDColumn(4, 0, x_pos);      
    man->AddNtupleRow(4);


     //hist->Draw();

    // Cleanup
    delete hist;
    delete gaussFit1;
    //hist->Draw();
}


G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory*)        
{

    SensorHit* aSensorHit = new SensorHit();

    G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();    
    G4Track *track = aStep->GetTrack();

    //track->SetTrackStatus(fStopAndKill);
    G4ParticleDefinition*  particle =aStep->GetTrack()->GetDefinition();
    if (particle->GetParticleName() == "opticalphoton"){   
       G4StepPoint *preStepPoint = aStep->GetPreStepPoint();//used when the neutron enters the detector
       G4StepPoint *postStepPoint = aStep->GetPostStepPoint();//used when the neutron leaves the detector
       const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();

       G4ThreeVector posPhotons = postStepPoint->GetPosition();//accessing the position       
       G4String particleName = particle->GetParticleName();
       //G4cout << "This is the particle in the sensor: " << particleName << G4endl;       

  
       G4int copyNo = touchable->GetCopyNumber();
       G4VPhysicalVolume *physVol = touchable->GetVolume();
       G4ThreeVector posDetector = physVol->GetTranslation();

       G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
       // Declare a vector to store the copy numbers
       std::vector<G4int> copyNumbers;

       G4String volumeName = touchable->GetVolume()->GetName(); // name of the volume
      // G4cout << "Photon detected in " << volumeName << "with Copy Number" << copyNo<< G4endl;
       G4AnalysisManager *man = G4AnalysisManager::Instance();

       G4double X1 = posPhotons.getX();
       G4double Y1 = posPhotons.getY();
       //std::vector<int> sensor1CopyNumbers;
      
        
       if (volumeName == "sensor_Vol1") {
          RecordSensorData(volumeName, 0, posDetector[0], posDetector[1], evt, copyNo, man);
          
          aSensorHit->SetSensorPosition(postStepPoint->GetPosition());
          aSensorHit->SetSensorEnergy(track->GetKineticEnergy());
          aSensorHit->SetSensorNumber(copyNo);
          if (SensorCollection) {
             HitID = SensorCollection->insert(aSensorHit);
          } else {
            G4cerr << "SensorCollection not initialized!" << G4endl;
            delete aSensorHit; // Avoid memory leaks
          }

          //} else {
            //delete aSensorHit; // Prevent memory leaks
          //}
           // G4cout << "Deposited energy: " << fTotalEnergyDeposited << G4endl;

          /*
          } else if (volumeName == "sensor_Vol2") {
            RecordSensorData(volumeName, 1, posDetector[0], posDetector[1], evt, copyNo, man);
            } else if (volumeName == "sensor_Vol3") {
              RecordSensorData(volumeName, 2, posDetector[0], posDetector[1], evt, copyNo, man);
              } else if (volumeName == "sensor_Vol4") {
                RecordSensorData(volumeName, 3, posDetector[0], posDetector[1], evt, copyNo, man);
                SaveToCSV(posDetector, copyNo);
                }*/

         //  FitHistogram(sensor1CopyNumbers);


       }
  }     


    return true;

}



void MySensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE) {
    // Log the total number of hits for debugging
    G4int nHits = SensorCollection->entries();
    G4cout << "End of Event: Number of hits in SensorCollection: " << nHits << G4endl;

    // Optionally process hits
    for (G4int i = 0; i < nHits; i++) {
        SensorHit* hit = (*SensorCollection)[i];
        /*G4cout << "Hit " << i << ": "
               << "Position = " << hit->GetSensorPosition()
               << ", Energy = " << hit->GetSensorEnergy()
               << ", Sensor ID = " << hit->GetSensorNumber() << G4endl;*/

         sensor1CopyNumbers.push_back(hit->GetSensorNumber());
    }

     // Fit the histogram with the collected data
    FitHistogram(sensor1CopyNumbers);
}


