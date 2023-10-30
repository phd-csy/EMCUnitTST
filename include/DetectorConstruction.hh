#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "CreateMapFromCSV.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "tls.hh"

#include <utility>
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4GlobalMagFieldMessenger;

class DetectorMessenger;

/// Detector construction class to define materials, geometry

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    DetectorConstruction();
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;
    auto GetCellNumber() const -> const auto& { return cellNumber; }
    // Set methods

private:
    // method

    // static data members

    // data members
    G4bool fCheckOverlaps = true; // option to activate checking of volumes overlaps
    G4int cellNumber;
};

#endif
