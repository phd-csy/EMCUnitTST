#include "SiPMSD.hh"
#include "SiPMHit.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"

SiPMSD::SiPMSD(const G4String &sdname, const G4String &hcName) : G4VSensitiveDetector(sdname) { collectionName.insert(hcName); }

SiPMSD::~SiPMSD() {}

void SiPMSD::Initialize(G4HCofThisEvent *hcOfThisEvent)
{
    G4int cellNumberTotal = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction())->GetCellNumber();

    hc = new SiPMHC(SensitiveDetectorName, collectionName[0]);
    auto hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    hcOfThisEvent->AddHitsCollection(hcID, hc);
}

G4bool SiPMSD::ProcessHits(G4Step *step, G4TouchableHistory *)
{
    auto particleDefinition = step->GetTrack()->GetDefinition();
    if (particleDefinition != G4OpticalPhoton::OpticalPhotonDefinition())
        return false;

    auto hit = new SiPMHit;
    auto globalTime = step->GetPostStepPoint()->GetGlobalTime();
    auto copyNo = step->GetTrack()->GetVolume()->GetCopyNo();
    hit->SetGlobalTime(globalTime);
    hit->SetCopyNo(copyNo);
    hc->insert(hit);

    return true;
}

void SiPMSD::EndOfEvent(G4HCofThisEvent*)
{
    // auto nofHits = fHitsCollection->entries();
}