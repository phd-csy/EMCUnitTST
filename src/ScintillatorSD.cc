#include "ScintillatorSD.hh"
#include "ScintillatorHit.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"

ScintillatorSD::ScintillatorSD(const G4String &sdName, const G4String &hcName)
: G4VSensitiveDetector(sdName)
{
    collectionName.insert(hcName);
}

ScintillatorSD::~ScintillatorSD()
{}

void ScintillatorSD::Initialize(G4HCofThisEvent *hcOfThisEvent)
{
    G4int cellNumberTotal = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction())->GetCellNumber();

    hc = new ScintillatorHC(SensitiveDetectorName, collectionName[0]);
    auto hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    hcOfThisEvent->AddHitsCollection(hcID, hc);
    for (auto i = 0; i < cellNumberTotal; i++)
        hc->insert(new ScintillatorHit());
}

G4bool ScintillatorSD::ProcessHits(G4Step *step, G4TouchableHistory *)
{
    auto copyNo = step->GetTrack()->GetVolume()->GetCopyNo();
    (*hc)[copyNo]->AddEnergyDeposit(step->GetTotalEnergyDeposit());

    return true;
}