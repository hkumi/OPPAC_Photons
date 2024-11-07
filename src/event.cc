#include "event.hh"

MyEventAction::MyEventAction(RunAction*)
{
    fEdep = 0.;
}

MyEventAction::~MyEventAction()
{}

void MyEventAction::BeginOfEventAction(const G4Event*)
{
    fEdep = 0.;
}

void MyEventAction::EndOfEventAction(const G4Event*)
{
    //G4cout << "Energy deposition: " << std::setw(7) << G4BestUnit(fEdep,"Energy") << G4endl;

    G4AnalysisManager *man = G4AnalysisManager::Instance();

    man->FillNtupleDColumn(4, 0, fEdep /CLHEP::keV);

    man->AddNtupleRow(4);
}
