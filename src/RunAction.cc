#include "RunAction.hh"

#include "DetectorConstruction.hh"
#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Timer.hh"
#include "PhotonHit.hh"
#include "PhotonSD.hh"
#include "ScintillatorHit.hh"
#include "ScintillatorSD.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction() :
    G4UserRunAction() {
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("root");

#ifdef G4MULTITHREADED
    analysisManager->SetNtupleMerging(true);
#endif

    // analysisManager->OpenFile("test");

    analysisManager->CreateNtuple("cellHit", "result");
    analysisManager->CreateNtupleIColumn("cellID");
    analysisManager->CreateNtupleDColumn("energyDeposit");
    analysisManager->CreateNtupleIColumn("nPhotons");
    analysisManager->FinishNtuple();

    analysisManager->CreateNtuple("pulse", "result");
    analysisManager->CreateNtupleIColumn("eventID");
    analysisManager->CreateNtupleIColumn("sipmID");
    analysisManager->CreateNtupleDColumn("nToF");
    analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* run) {
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*) {
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
