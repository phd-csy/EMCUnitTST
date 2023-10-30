#ifndef SIPMHIT_HH
#define SIPMHIT_HH

#include "G4THitsCollection.hh"
#include "G4VHit.hh"

class SiPMHit : public G4VHit {
public:
    void SetGlobalTime(G4double time) { globalTime = time; }
    G4double GetGlobalTime() const { return globalTime; }

    void SetCopyNo(G4int copyNumber) { copyNo = copyNumber; }
    G4double GetCopyNo() const { return copyNo; }

private:
    G4double globalTime;
    G4int copyNo;
};

using SiPMHC = G4THitsCollection<SiPMHit>;

#endif