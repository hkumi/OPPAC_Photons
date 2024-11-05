#ifndef StepAction_h
#define StepAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <vector>
#include "DC.hh"  // Include the DetectorConstruction header
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
class StepAction : public G4UserSteppingAction
{

public:
  StepAction(G4String data_file, double lunghezza_collimatore);
  ~StepAction();

  virtual  void UserSteppingAction(const G4Step*);
           double hysto(double valore);
  
private:
  G4String datai_file;
  //std::vector<G4double> coor_x1;
  //std::vector<G4double> coor_x2; 
  //std::vector<G4double> coor_y1; 
  //std::vector<G4double> coor_y2;  
  int s_x1[33];
  int s_x2[33];
  int s_y1[33];
  int s_y2[33];
  double collimatore;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
