#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>  
#include "G4AnalysisManager.hh"
#include "SensorHit.hh"
#include "G4VHitsCollection.hh"
#include "G4EventManager.hh"
class G4HCofThisEvent;
class MySensitiveDetector : public G4VSensitiveDetector
{
public:
    MySensitiveDetector(G4String);
    ~MySensitiveDetector();

    double calculateWeightedMeanX(double Px1, double Nx1, double sigmax1);    
    void SaveToCSV(G4ThreeVector posPhotons, G4int copyNo);
    void RecordSensorData(const std::string& volumeName, int ntupleIndex, double posX, double posY, int event, int copyNo, G4AnalysisManager* man);  
    void FitHistogram(const std::vector<double>& Position);
    std::vector<int>& GetSensor1CopyNumbers()  { return sensor1CopyNumbers; }
    G4bool ProcessHits(G4Step *, G4TouchableHistory*);
    void Initialize(G4HCofThisEvent*);
    void EndOfEvent(G4HCofThisEvent*);

    
private:
    //virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);
    std::vector<int> sensor1CopyNumbers;
    
    //virtual void Initialize(G4HCofThisEvent*) override;
    //virtual void EndOfEvent(G4HCofThisEvent *) override;

    G4double fTotalEnergyDeposited;
    G4int HitID;
    SensorHitsCollection* SensorCollection;
    

};

#endif

