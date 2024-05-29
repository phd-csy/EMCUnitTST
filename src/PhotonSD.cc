#include "PhotonSD.hh"

#include "DetectorConstruction.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4RunManager.hh"
#include "PhotonHit.hh"

PhotonSD::PhotonSD(const G4String& sdname, const G4String& hcName) :
    G4VSensitiveDetector(sdname) { collectionName.insert(hcName); }

PhotonSD::~PhotonSD() {}

void PhotonSD::Initialize(G4HCofThisEvent* hcOfThisEvent) {
    G4int cellNumberTotal = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction())->GetCellNumber();

    hc = new PhotonHC(SensitiveDetectorName, collectionName[0]);
    auto hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    hcOfThisEvent->AddHitsCollection(hcID, hc);
}

G4bool PhotonSD::ProcessHits(G4Step* step, G4TouchableHistory*) {
    auto particleDefinition = step->GetTrack()->GetDefinition();
    if (particleDefinition == G4OpticalPhoton::Definition()) {

        step->GetTrack()->SetTrackStatus(fStopAndKill);

        auto hit = new PhotonHit;
        auto globalTime = step->GetPostStepPoint()->GetGlobalTime();
        auto copyNo = step->GetTrack()->GetVolume()->GetCopyNo();
        hit->SetGlobalTime(globalTime);
        hit->SetCopyNo(copyNo);
        hc->insert(hit);

        return true;
    }
    return false;
}

void PhotonSD::EndOfEvent(G4HCofThisEvent*) {}