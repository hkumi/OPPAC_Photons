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
 //   G4cout << "SensitiveDetectorName: " << SensitiveDetectorName << G4endl;
//G4cout << "Collection Name: " << collectionName[0] << G4endl;

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
double MySensitiveDetector::calculateWeightedMeanX(double P_x1, double N_x1, double sigma_x1, 
                                                   double P_x2, double N_x2, double sigma_x2) {
    // Check for invalid sigma values
    if (sigma_x1 == 0 || sigma_x2 == 0) {
        std::cout << "Error: Sigma values cannot be zero!" << std::endl;
        return -1; // Return a sentinel value to indicate error
    }

    // Calculate the weighted mean
    double numerator = (P_x1 * N_x1 / sigma_x1) + (P_x2 * N_x2 / sigma_x2);
    double denominator = (N_x1 / sigma_x1) + (N_x2 / sigma_x2);

    return numerator / denominator;
}

// Calculate the weighted mean for the y-coordinate
double MySensitiveDetector::calculateWeightedMeanY(double P_y1, double N_y1, double sigma_y1, 
                                                   double P_y2, double N_y2, double sigma_y2) {
    // Check for invalid sigma values
    if (sigma_y1 == 0 || sigma_y2 == 0) {
        std::cout << "Error: Sigma values cannot be zero!" << std::endl;
        return -1; // Return a sentinel value to indicate error
    }

    // Calculate the weighted mean
    double numerator = (P_y1 * N_y1 / sigma_y1) + (P_y2 * N_y2 / sigma_y2);
    double denominator = (N_y1 / sigma_y1) + (N_y2 / sigma_y2);

    return numerator / denominator;
}




void MySensitiveDetector::FitHistogram(const std::vector<double>& Position1, const std::vector<double>& Position2) {
    if (Position1.empty() || Position2.empty()) return;

    // Find min and max values in Position1 and Position2 vectors
    double min1 = *std::min_element(Position1.begin(), Position1.end());
    double max1 = *std::max_element(Position1.begin(), Position1.end());
    double min2 = *std::min_element(Position2.begin(), Position2.end());
    double max2 = *std::max_element(Position2.begin(), Position2.end());

    int bins = 50;

    // Create histograms with unique names
    static int histCounter1 = 0, histCounter2 = 0;
    TH1D* hist1 = new TH1D(Form("copy1_%d", histCounter1++), "Sensor1 Distribution", bins, min1, max1);
    TH1D* hist2 = new TH1D(Form("copy2_%d", histCounter2++), "Sensor2 Distribution", bins, min2, max2);

    for (double Pos1 : Position1) hist1->Fill(Pos1);
    for (double Pos2 : Position2) hist2->Fill(Pos2);

    // Fit histograms with Gaussian functions
    TF1* gaussFit1 = new TF1("gaussFit1", "gaus", min1, max1);
    gaussFit1->SetLineColor(kRed);
    hist1->Fit(gaussFit1, "R");

    TF1* gaussFit2 = new TF1("gaussFit2", "gaus", min2, max2);
    gaussFit2->SetLineColor(kRed);
    hist2->Fit(gaussFit2, "R");

    // Extract Gaussian fit parameters
    double mean1 = gaussFit1->GetParameter(1);
    double sigma1 = gaussFit1->GetParameter(2);
    double amplitude1 = gaussFit1->GetParameter(0);

    double mean2 = gaussFit2->GetParameter(1);
    double sigma2 = gaussFit2->GetParameter(2);
    double amplitude2 = gaussFit2->GetParameter(0);

    // Debugging outputs
    G4cout << "Sensor1 Gaussian Fit: Mean=" << mean1 << ", Sigma=" << sigma1
           << ", Amplitude=" << amplitude1 << G4endl;
    G4cout << "Sensor2 Gaussian Fit: Mean=" << mean2 << ", Sigma=" << sigma2
           << ", Amplitude=" << amplitude2 << G4endl;

    // Calculate weighted mean
    double y_pos = calculateWeightedMeanY(mean1, amplitude1, sigma1, mean2, amplitude2, sigma2);
    G4cout << "Weighted Mean X: " << y_pos << G4endl;

    // Store the result in ntuple
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    man->FillNtupleDColumn(4, 0, y_pos/mm);
    man->AddNtupleRow(4);

    // Cleanup
    delete hist1;
    delete hist2;
    delete gaussFit1;
    delete gaussFit2;
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
          aSensorHit->SetVolumeName(volumeName);

          if (SensorCollection) {
             HitID = SensorCollection->insert(aSensorHit);
            
          } else {
            G4cerr << "SensorCollection not initialized!" << G4endl;
            delete aSensorHit; // Avoid memory leaks
            }
       } else if (volumeName == "sensor_Vol2") {
         // Create a new hit object for sensor_Vol2
         SensorHit* anotherSensorHit = new SensorHit();

         RecordSensorData(volumeName, 1, posDetector[0], posDetector[1], evt, copyNo, man);

         anotherSensorHit->SetSensorPosition(postStepPoint->GetPosition());
         anotherSensorHit->SetSensorEnergy(track->GetKineticEnergy());
         anotherSensorHit->SetVolumeName(volumeName);

         if (SensorCollection) {
            HitID = SensorCollection->insert(anotherSensorHit);

         } else {
           G4cerr << "SensorCollection not initialized!" << G4endl;
           delete anotherSensorHit; // Avoid memory leaks
           }
       }

         /*
            } else if (volumeName == "sensor_Vol3") {
              RecordSensorData(volumeName, 2, posDetector[0], posDetector[1], evt, copyNo, man);
              } else if (volumeName == "sensor_Vol4") {
                RecordSensorData(volumeName, 3, posDetector[0], posDetector[1], evt, copyNo, man);
                SaveToCSV(posDetector, copyNo);
                }*/

         //  FitHistogram(sensor1CopyNumbers);


       
    }     


    return true;

}


void MySensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE) {
    // Log the total number of hits for debugging
    G4int nHits = SensorCollection->entries();
    //G4cout << "End of Event: Number of hits in SensorCollection: " << nHits << G4endl;

    // Check if there are any hits
    if (nHits > 0) {
        std::vector<double> yPositions1, yPositions2; // For sensors 1 and 2
        std::vector<double> xPositions1, xPositions2;
        std::vector<G4int> sensor1CopyNumbers, sensor2CopyNumbers;

        // Process hits
        for (G4int i = 0; i < nHits; i++) {
            SensorHit* hit = (*SensorCollection)[i];
            G4ThreeVector position = hit->GetSensorPosition(); // Get the hit position
            G4String VolumeName = hit->GetVolumeName();       // Copy number (optional)

            // Separate data based on the sensor volume name
            if (VolumeName == "sensor_Vol1") { // For sensor_Vol1
                xPositions1.push_back(position.x());
                yPositions1.push_back(position.y());
                //sensor1CopyNumbers.push_back(sensorNumber);
            } else if (VolumeName == "sensor_Vol2") { // For sensor_Vol2
                xPositions2.push_back(position.x());
                yPositions2.push_back(position.y());
                //sensor2CopyNumbers.push_back(sensorNumber);
            }
        }

        // Output for debugging
        G4cout << "Sensor_Vol1 Hits: " << yPositions1.size() 
               << ", Sensor_Vol2 Hits: " << yPositions2.size() << G4endl;

        // Perform Gaussian fits or any other processing
        if (!yPositions1.empty() ||!yPositions2.empty()) {
            FitHistogram(yPositions1,yPositions2); // Fit data for Sensor_Vol1 and Sensor_Vol2
        }
        
    }
}

