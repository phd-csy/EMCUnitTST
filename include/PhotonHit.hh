#ifndef PHOTONHIT_HH
#define PHOTONHIT_HH

#include "G4THitsCollection.hh"
#include "G4VHit.hh"

class PhotonHit : public G4VHit {
public:
    void SetGlobalTime(G4double time) { globalTime = time; }
    G4double GetGlobalTime() const { return globalTime; }

    void SetCopyNo(G4int copyNumber) { copyNo = copyNumber; }
    G4double GetCopyNo() const { return copyNo; }

private:
    G4double globalTime;
    G4int copyNo;
};

using PhotonHC = G4THitsCollection<PhotonHit>;

#endif