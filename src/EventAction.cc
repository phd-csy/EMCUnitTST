#include "EventAction.hh"

#include "DetectorConstruction.hh"
#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4Trajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "PhotonHit.hh"
#include "PhotonSD.hh"
#include "ScintillatorHit.hh"
#include "ScintillatorSD.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event) {
    auto scintillatorHCid = G4SDManager::GetSDMpointer()->GetCollectionID("ScintillatorHitsCollection");
    auto scintHC = static_cast<ScintillatorHC*>(event->GetHCofThisEvent()->GetHC(scintillatorHCid));

    auto photonHCid = G4SDManager::GetSDMpointer()->GetCollectionID("PhotonHitsCollection");
    auto photonHC = static_cast<PhotonHC*>(event->GetHCofThisEvent()->GetHC(photonHCid));

    auto analysisManager = G4AnalysisManager::Instance();

    G4int cellNumberTotal = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction())->GetCellNumber();
    auto eventID = event->GetEventID();

    for (auto i = 0; i < cellNumberTotal; i++) {

        auto energyDeposit = (*scintHC)[i]->GetEnergyDeposit();

        if (energyDeposit > 0.) {
            analysisManager->FillNtupleIColumn(0, 0, i);
            analysisManager->FillNtupleDColumn(0, 1, energyDeposit);
            analysisManager->FillNtupleIColumn(0, 2, photonHC->entries());
            analysisManager->AddNtupleRow(0);
        }
    }

    for (long unsigned int j = 0; j < photonHC->entries(); j++) {

        auto photonHit = (*photonHC)[j];

        analysisManager->FillNtupleIColumn(1, 0, eventID);
        analysisManager->FillNtupleIColumn(1, 1, photonHit->GetCopyNo());
        analysisManager->FillNtupleDColumn(1, 2, photonHit->GetGlobalTime());
        analysisManager->AddNtupleRow(1);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
