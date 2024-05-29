#ifndef PHOTONSD_HH
#define PHOTONSD_HH

#include "G4OpticalPhoton.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "PhotonHit.hh"

class PhotonSD : public G4VSensitiveDetector {
public:
    PhotonSD(const G4String&, const G4String&);
    ~PhotonSD() override;

    void Initialize(G4HCofThisEvent*) override;
    G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;
    void EndOfEvent(G4HCofThisEvent*) override;

private:
    PhotonHC* hc;
};

#endif