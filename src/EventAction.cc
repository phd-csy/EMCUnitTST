#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "ScintillatorHit.hh"
#include "ScintillatorSD.hh"
#include "SiPMHit.hh"
#include "SiPMSD.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
    auto scintillatorHCid = G4SDManager::GetSDMpointer()->GetCollectionID("ScintillatorHitsCollection");
    auto scintHC = static_cast<ScintillatorHC *>(event->GetHCofThisEvent()->GetHC(scintillatorHCid));

    auto sipmHCid = G4SDManager::GetSDMpointer()->GetCollectionID("SiPMHitsCollection");
    auto sipmHC = static_cast<SiPMHC *>(event->GetHCofThisEvent()->GetHC(sipmHCid));

    auto analysisManager = G4AnalysisManager::Instance();

    G4int cellNumberTotal = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction())->GetCellNumber();
    auto eventID = event->GetEventID();

    for (auto i = 0; i < cellNumberTotal; i++){

        auto energyDeposit = (*scintHC)[i]->GetEnergyDeposit();

        if (energyDeposit > 0.){
            analysisManager->FillNtupleIColumn(0, 0, i);
            analysisManager->FillNtupleDColumn(0, 1, energyDeposit);
            analysisManager->FillNtupleIColumn(0, 2, sipmHC->GetSize());
            analysisManager->AddNtupleRow(0);
        }
    }

    for (long unsigned int j = 0; j < sipmHC->GetSize(); j++){

        auto sipmHit = (*sipmHC)[j];

        analysisManager->FillNtupleIColumn(1, 0, eventID);
        analysisManager->FillNtupleIColumn(1, 1, sipmHit->GetCopyNo());
        analysisManager->FillNtupleDColumn(1, 2, sipmHit->GetGlobalTime());
        analysisManager->AddNtupleRow(1);
    }

    auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
    if ((printModulo > 0) && (eventID % printModulo == 0))
    {
        G4cout << "---> End of event: " << eventID << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


