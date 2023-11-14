#include "DetectorConstruction.hh"

#include "G4AutoDelete.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Element.hh"
#include "G4ExtrudedSolid.hh"
#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4Isotope.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4OpticalMaterialProperties.hh"
#include "G4OpticalSurface.hh"
#include "G4PhysicalConstants.hh"
#include "G4PVPlacement.hh"
#include "G4QuadrangularFacet.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4TessellatedSolid.hh"
#include "G4ThreeVector.hh"
#include "G4Trd.hh"
#include "G4TriangularFacet.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "pmp/algorithms/differential_geometry.h"
#include "pmp/algorithms/normals.h"
#include "pmp/algorithms/subdivision.h"
#include "pmp/algorithms/utilities.h"
#include "pmp/surface_mesh.h"
#include "ScintillatorHit.hh"
#include "ScintillatorSD.hh"
#include "SiPMHit.hh"
#include "SiPMSD.hh"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <utility>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::DetectorConstruction() {
    cellNumber = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct() {

    // Define Material

    const auto nistManager = G4NistManager::Instance();

    const auto galactic = new G4Material("galactic", 1, 1.008 * g / mole, 1.e-25 * g / cm3, kStateGas, 2.73 * kelvin, 3.e-18 * pascal);
    const auto air = nistManager->FindOrBuildMaterial("G4_AIR");
    const auto aluminum = nistManager->FindOrBuildMaterial("G4_Al");
    const auto silicon = nistManager->FindOrBuildMaterial("G4_Si");
    const auto pvc = nistManager->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");

    const auto carbonElement = nistManager->FindOrBuildElement("C");
    const auto hydrogenElement = nistManager->FindOrBuildElement("H");
    const auto oxygenElement = nistManager->FindOrBuildElement("O");
    const auto siliconElement = nistManager->FindOrBuildElement("Si");
    const auto cesiumElement = nistManager->FindOrBuildElement("Cs");
    const auto iodideElement = nistManager->FindOrBuildElement("I");
    const auto thaliumElement = nistManager->FindOrBuildElement("Tl");
    const auto bromideElement = nistManager->FindOrBuildElement("Br");
    const auto lanthanumElement = nistManager->FindOrBuildElement("La");
    const auto yttriumElement = nistManager->FindOrBuildElement("Y");
    const auto lutetiumElement = nistManager->FindOrBuildElement("Lu");

    const auto siliconeOil = new G4Material("silicone_oil", 0.97 * g / cm3, 4, kStateLiquid);
    siliconeOil->AddElement(carbonElement, 2);
    siliconeOil->AddElement(hydrogenElement, 6);
    siliconeOil->AddElement(oxygenElement, 1);
    siliconeOil->AddElement(siliconElement, 1);

    const auto glass = new G4Material("Fused Silica", 2.64 * g / cm3, 2, kStateSolid);
    glass->AddElement(oxygenElement, 0.532570);
    glass->AddElement(siliconElement, 0.467430);

    const auto csI = new G4Material("CsI", 4.51 * g / cm3, 3, kStateSolid);
    csI->AddElement(cesiumElement, 0.507556);
    csI->AddElement(iodideElement, 0.484639);
    csI->AddElement(thaliumElement, 0.007805);

    // const auto labr = new G4Material("LaBr3", 5.08 * g / cm3, 3, kStateSolid);
    // labr->AddElement(bromideElement, 0.631308);
    // labr->AddElement(cesiumElement, 0.036582);
    // labr->AddElement(lanthanumElement, 0.332110);

    // const auto lyso = new G4Material("LYSO", 7.1 * g / cm3, 5, kStateSolid);
    // lyso->AddElement(oxygenElement, 0.175801);
    // lyso->AddElement(siliconElement, 0.061720);
    // lyso->AddElement(yttriumElement, 0.019538);
    // lyso->AddElement(lutetiumElement, 0.730562);
    // lyso->AddElement(cesiumElement, 0.012379);

    /////////////////////////////////////////////
    // Construct Material Properties Tables
    /////////////////////////////////////////////

    const auto fLambda_min = 200 * nm;
    const auto fLambda_max = 700 * nm;
    std::vector<G4double> ENERGY_BIN = {h_Planck * c_light / fLambda_max, h_Planck * c_light / fLambda_min};

    //============================================ air ================================================

    std::vector<G4double> RINDEX_AIR = {1.00, 1.00};
    // Air refractive index at 20 oC and 1 atm (from PDG)
    for (auto&& i : RINDEX_AIR) {
        i = i + 2.73 * std::pow(10.0, -4);
    }

    auto MPT_Air = new G4MaterialPropertiesTable();
    MPT_Air->AddProperty("RINDEX", ENERGY_BIN, RINDEX_AIR);
    galactic->SetMaterialPropertiesTable(MPT_Air);

    //============================================ silicone oil =======================================

    // silicone oil refractive index
    std::vector<G4double> ENERGY_SiOIL = {8.01532E-07, 1.89386E-06, 1.92915E-06, 2.1093E-06, 2.27541E-06, 2.55633E-06, 2.58828E-06};
    std::vector<G4double> RINDEX_SiOIL = {1.3912, 1.3992, 1.3997, 1.4015, 1.4034, 1.4071, 1.4076};

    // silicone oil absoption length
    std::vector<G4double> ABSL_SiOIL = {15 * cm, 15 * cm};

    auto siliconeOilPropertiesTable = new G4MaterialPropertiesTable();
    siliconeOilPropertiesTable->AddProperty("RINDEX", ENERGY_SiOIL, RINDEX_SiOIL);
    siliconeOilPropertiesTable->AddProperty("ABSLENGTH", ENERGY_BIN, ABSL_SiOIL);
    siliconeOil->SetMaterialPropertiesTable(siliconeOilPropertiesTable);

    //============================================ Quartz =============================================

    std::vector<G4double> RINDEX_QUARTZ = {1.54, 1.54};

    auto windowPropertiesTable = new G4MaterialPropertiesTable();
    windowPropertiesTable->AddProperty("RINDEX", ENERGY_BIN, RINDEX_QUARTZ);
    glass->SetMaterialPropertiesTable(windowPropertiesTable);
    windowPropertiesTable->DumpTable();

    //============================================ Crystal ============================================

    // CsI(Tl)
    std::vector<G4double> GROUPV_CSI = {167.482, 167.482};

    auto csiPropertiesTable = new G4MaterialPropertiesTable;
    auto csiProperties(CreateMapFromCSV<G4double>("data/CsI_properties.csv"));
    csiPropertiesTable->AddProperty(
        "RINDEX",
        csiProperties["energy"],
        csiProperties["RINDEX"]);
    csiPropertiesTable->AddProperty(
        "GROUPVEL",
        ENERGY_BIN,
        GROUPV_CSI);
    csiPropertiesTable->AddProperty(
        "ABSLENGTH",
        csiProperties["energy"],
        csiProperties["ABSLENGTH"]);
    csiPropertiesTable->AddProperty(
        "SCINTILLATIONCOMPONENT1",
        csiProperties["energy"],
        csiProperties["SCINTILLATIONCOMPONENT1"]);
    csiPropertiesTable->AddConstProperty("SCINTILLATIONYIELD", csiProperties["SCINTILLATIONYIELD"].front());
    csiPropertiesTable->AddConstProperty("SCINTILLATIONTIMECONSTANT1", csiProperties["SCINTILLATIONTIMECONSTANT1"].front());
    csiPropertiesTable->AddConstProperty("RESOLUTIONSCALE", csiProperties["RESOLUTIONSCALE"].front());
    csI->SetMaterialPropertiesTable(csiPropertiesTable);
    csiPropertiesTable->DumpTable();

    // // LaBr3(Ce)
    // auto labrPropertiesTable = new G4MaterialPropertiesTable;
    // auto labrProperties(CreateMapFromCSV<G4double>("data/LaBr_properties.csv"));
    // labrPropertiesTable->AddProperty(
    //     "RINDEX",
    //     labrProperties["energy"],
    //     labrProperties["RINDEX"]);
    // labrPropertiesTable->AddProperty(
    //     "ABSLENGTH",
    //     labrProperties["energy"],
    //     labrProperties["ABSLENGTH"]);
    // labrPropertiesTable->AddProperty(
    //     "SCINTILLATIONCOMPONENT1",
    //     labrProperties["energy"],
    //     labrProperties["SCINTILLATIONCOMPONENT1"]);
    // labrPropertiesTable->AddConstProperty("SCINTILLATIONYIELD", labrProperties["SCINTILLATIONYIELD"].front());
    // labrPropertiesTable->AddConstProperty("SCINTILLATIONTIMECONSTANT1", labrProperties["SCINTILLATIONTIMECONSTANT1"].front());
    // labrPropertiesTable->AddConstProperty("RESOLUTIONSCALE", labrProperties["RESOLUTIONSCALE"].front());
    // labr->SetMaterialPropertiesTable(labrPropertiesTable);

    // // LYSO(Ce)
    // auto lysoPropertiesTable = new G4MaterialPropertiesTable;
    // auto lysoProperties(CreateMapFromCSV<G4double>("data/LYSO_properties.csv"));
    // lysoPropertiesTable->AddProperty(
    //     "RINDEX",
    //     lysoProperties["energy"],
    //     lysoProperties["RINDEX"]);
    // lysoPropertiesTable->AddProperty(
    //     "ABSLENGTH",
    //     lysoProperties["energy"],
    //     lysoProperties["ABSLENGTH"]);
    // lysoPropertiesTable->AddProperty(
    //     "SCINTILLATIONCOMPONENT1",
    //     lysoProperties["energy"],
    //     lysoProperties["SCINTILLATIONCOMPONENT1"]);
    // lysoPropertiesTable->AddConstProperty("SCINTILLATIONYIELD", lysoProperties["SCINTILLATIONYIELD"].front());
    // lysoPropertiesTable->AddConstProperty("SCINTILLATIONTIMECONSTANT1", lysoProperties["SCINTILLATIONTIMECONSTANT1"].front());
    // lysoPropertiesTable->AddConstProperty("RESOLUTIONSCALE", lysoProperties["RESOLUTIONSCALE"].front());
    // lyso->SetMaterialPropertiesTable(lysoPropertiesTable);

    // Define World Volume

    const auto worldSize = 2 * m;
    const auto worldS = new G4Box("world", 0.5 * worldSize, 0.5 * worldSize, 0.5 * worldSize);
    const auto worldLV = new G4LogicalVolume(worldS, galactic, "World");
    const auto worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World", nullptr, false, 0, true);

    const auto fCrystalHalfWidth = 30 * mm;
    const auto fCrystalHalfLength = 7.5 * cm;

    const auto fCouplerHalfWidth = 35 * mm;
    const auto fCouplerHalfThickness = 200 * um;

    const auto fPMTWindowDiameter = fCouplerHalfWidth;
    const auto fPMTWindowHalfThickness = 0.5 * mm;

    const auto fSiPMHalfWidth = fCouplerHalfWidth;
    const auto fSiPMHalfThickness = 10 * nm;

    const auto offsetlength = 0 * cm;

    const auto translation1 = G4ThreeVector(0., offsetlength, fCrystalHalfLength + fCouplerHalfThickness);
    const auto translation2 = G4ThreeVector(0., offsetlength, fCrystalHalfLength + fCouplerHalfThickness * 2 + fPMTWindowHalfThickness);
    const auto translation3 = G4ThreeVector(0., offsetlength, fCrystalHalfLength + fCouplerHalfThickness * 2 + fPMTWindowHalfThickness * 2 + fSiPMHalfThickness);

    // Construct Volume

    const G4int nsect = 5;
    std::vector<G4TwoVector> polygon(nsect);
    G4double ang = twopi / nsect;
    for (G4int i = 0; i < nsect; ++i) {
        G4double phi = i * ang;
        G4double cosphi = std::cos(phi);
        G4double sinphi = std::sin(phi);
        polygon[i].set(fCrystalHalfWidth * cosphi, fCrystalHalfWidth * sinphi);
    }

    G4TwoVector offsetA(0, 0), offsetB(0, offsetlength);
    G4double scaleA = 0.5, scaleB = 1;

    const auto crystalSV = new G4ExtrudedSolid("Extruded", polygon, fCrystalHalfLength, offsetA, scaleA, offsetB, scaleB);

    // const auto crystalSV = new G4Box("crystal", fCrystalHalfWidth, fCrystalHalfWidth, fCrystalHalfLength);
    // const auto crystalSV = new G4Trd("crystal", fCrystalHalfWidth - 20 * mm, fCrystalHalfWidth, fCrystalHalfWidth - 20 * mm, fCrystalHalfWidth, fCrystalHalfLength);

    //========================================== CsI(Tl) ============================================
    const auto crystalLV = new G4LogicalVolume{crystalSV, csI, "crystal"};
    //========================================== LaBr3(Ce) ==========================================
    // const auto crystalLV = new G4LogicalVolume{crystalSV, labr, "crystal"};
    //========================================== LYSO(Ce) ===========================================
    // const auto crystalLV = new G4LogicalVolume{crystalSV, lyso, "crystal"};
    //===============================================================================================

    const auto crystalPV = new G4PVPlacement(nullptr, G4ThreeVector(), crystalLV, "crystal", worldLV, false, 0, true);

    const auto couplerSV = new G4Tubs("coupler", 0, fCouplerHalfWidth, fCouplerHalfThickness, 0, 2 * pi);
    const auto couplerLV = new G4LogicalVolume(couplerSV, siliconeOil, "coupler");
    const auto couplerPV = new G4PVPlacement(nullptr, translation1, couplerLV, "coupler", worldLV, false, 0, true);

    const auto windowSV = new G4Tubs("window", 0, fPMTWindowDiameter, fPMTWindowHalfThickness, 0, 2 * pi);
    const auto windowLV = new G4LogicalVolume(windowSV, glass, "window");
    const auto windowPV = new G4PVPlacement{nullptr, translation2, windowLV, "window", worldLV, false, 0, true};

    const auto sipmSV = new G4Tubs("sipm", 0, fSiPMHalfWidth, fSiPMHalfThickness, 0, 2 * pi);
    const auto sipmLV = new G4LogicalVolume(sipmSV, silicon, "sipm");
    new G4PVPlacement(nullptr, translation3, sipmLV, "sipm", worldLV, false, 0, true);

    // Define Surface

    std::vector<G4double> REFLECTIVITY = {0.955, 0.955};
    auto tempPropertiesTable = new G4MaterialPropertiesTable();
    tempPropertiesTable->AddProperty("REFLECTIVITY", ENERGY_BIN, REFLECTIVITY);

    auto alSurface = new G4OpticalSurface("Al foil", glisur, ground, dielectric_metal);
    new G4LogicalBorderSurface("AlSkinSurface", crystalPV, worldPV, alSurface);
    alSurface->SetMaterialPropertiesTable(tempPropertiesTable);

    auto optocouplerSurface = new G4OpticalSurface("Optocoupler", glisur, polished, dielectric_metal);
    new G4LogicalBorderSurface("optocouplerSurface", couplerPV, worldPV, optocouplerSurface);
    optocouplerSurface->SetMaterialPropertiesTable(tempPropertiesTable);

    auto windowSurface = new G4OpticalSurface("Optocoupler", glisur, polished, dielectric_metal);
    new G4LogicalBorderSurface("windowSurface", windowPV, worldPV, windowSurface);
    windowSurface->SetMaterialPropertiesTable(tempPropertiesTable);

    auto cathodeSurfacePropertiesTable = new G4MaterialPropertiesTable();
    auto cathodeSurfaceProperties(CreateMapFromCSV<G4double>("data/PMT_properties.csv"));

    std::vector<G4double> cathodeSurfacePropertiesEnergy(cathodeSurfaceProperties["wavelength"].size());
    std::vector<G4double> cathodeSurfacePropertiesEfficiency(cathodeSurfaceProperties["EFFICIENCY"].size());

    std::transform(cathodeSurfaceProperties["wavelength"].begin(), cathodeSurfaceProperties["wavelength"].end(), cathodeSurfacePropertiesEnergy.begin(),
                   [](auto val) { return h_Planck * c_light / (val * nm / mm); });

    std::transform(cathodeSurfaceProperties["EFFICIENCY"].begin(), cathodeSurfaceProperties["EFFICIENCY"].end(), cathodeSurfacePropertiesEfficiency.begin(),
                   [](auto n) { return n * perCent; });

    auto cathodeSurface = new G4OpticalSurface("Cathode", glisur, polished, dielectric_metal);
    new G4LogicalSkinSurface("sipmSkinSurface", sipmLV, cathodeSurface);
    cathodeSurfacePropertiesTable->AddProperty(
        "REFLECTIVITY",
        cathodeSurfacePropertiesEnergy,
        cathodeSurfaceProperties["REFLECTIVITY"]);
    cathodeSurfacePropertiesTable->AddProperty(
        "EFFICIENCY",
        cathodeSurfacePropertiesEnergy,
        cathodeSurfacePropertiesEfficiency);
    cathodeSurface->SetMaterialPropertiesTable(cathodeSurfacePropertiesTable);
    cathodeSurfacePropertiesTable->DumpTable();

    ++cellNumber;

    return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField() {
    // Sensitive detectors
    auto scintSD = new ScintillatorSD("ScintillatorSD", "ScintillatorHitsCollection");
    G4SDManager::GetSDMpointer()->AddNewDetector(scintSD);
    SetSensitiveDetector("crystal", scintSD, true);

    auto sipmSD = new SiPMSD("SiPMSD", "SiPMHitsCollection");
    G4SDManager::GetSDMpointer()->AddNewDetector(sipmSD);
    SetSensitiveDetector("sipm", sipmSD, true);
}