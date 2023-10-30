#ifndef SIPMSD_HH
#define SIPMSD_HH

#include "G4OpticalPhoton.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "SiPMHit.hh"

class SiPMSD : public G4VSensitiveDetector {
public:
    SiPMSD(const G4String&, const G4String&);
    ~SiPMSD() override;

    void Initialize(G4HCofThisEvent*) override;
    G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;
    void EndOfEvent(G4HCofThisEvent*) override;

private:
    SiPMHC* hc;
};

#endif