#ifndef SCINTILLATORHIT_HH
#define SCINTILLATORHIT_HH

#include "G4THitsCollection.hh"
#include "G4VHit.hh"

class ScintillatorHit : public G4VHit {
public:
    void AddEnergyDeposit(G4double eDep) { energyDeposit += eDep; }
    G4double GetEnergyDeposit() const { return energyDeposit; }

private:
    G4double energyDeposit = 0.;
};

using ScintillatorHC = G4THitsCollection<ScintillatorHit>;

#endif