#ifndef SCINTILLATORSD_HH
#define SCINTILLATORSD_HH

#include "G4VSensitiveDetector.hh"
#include "ScintillatorHit.hh"
#include "G4SDManager.hh"

class ScintillatorSD : public G4VSensitiveDetector
{
public:
    ScintillatorSD(const G4String &, const G4String &);
    ~ScintillatorSD() override;

    void Initialize(G4HCofThisEvent *) override;
    G4bool ProcessHits(G4Step *, G4TouchableHistory *) override;

private:
    ScintillatorHC *hc;
};

#endif