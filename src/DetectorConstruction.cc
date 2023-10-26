#include "DetectorConstruction.hh"
#include "ScintillatorHit.hh"
#include "ScintillatorSD.hh"
#include "SiPMHit.hh"
#include "SiPMSD.hh"

#include "G4AutoDelete.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4ExtrudedSolid.hh"
#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4PVPlacement.hh"
#include "G4QuadrangularFacet.hh"
#include "G4SDManager.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4TessellatedSolid.hh"
#include "G4ThreeVector.hh"
#include "G4TriangularFacet.hh"
#include "G4Tubs.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4RunManager.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"

#include "pmp/algorithms/differential_geometry.h"
#include "pmp/algorithms/normals.h"
#include "pmp/algorithms/subdivision.h"
#include "pmp/algorithms/utilities.h"
#include "pmp/surface_mesh.h"

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
    // const auto bgo = nistManager->FindOrBuildMaterial("G4_BGO");
    const auto silicon = nistManager->FindOrBuildMaterial("G4_Si");
    // const auto csI = nistManager->FindOrBuildMaterial("G4_CESIUM_IODIDE");
    const auto pvc = nistManager->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");

    const auto carbonElement = nistManager->FindOrBuildElement("C");
    const auto hydrogenElement = nistManager->FindOrBuildElement("H");
    const auto oxygenElement = nistManager->FindOrBuildElement("O");
    const auto siliconElement = nistManager->FindOrBuildElement("Si");

    const auto siliconeOil = new G4Material("silicone_oil", 0.97 * g / cm3, 4, kStateLiquid);
    siliconeOil->AddElement(carbonElement, 2);
    siliconeOil->AddElement(hydrogenElement, 6);
    siliconeOil->AddElement(oxygenElement, 1);
    siliconeOil->AddElement(siliconElement, 1);

    const auto cesiumElement = nistManager->FindOrBuildElement("Ce");
    const auto iodideElement = nistManager->FindOrBuildElement("I");
    const auto thaliumElement = nistManager->FindOrBuildElement("Tl");

    const auto bromideElement = nistManager->FindOrBuildElement("Br");
    const auto lanthanumElement = nistManager->FindOrBuildElement("La");

    const auto yttriumElement = nistManager->FindOrBuildElement("Y");
    const auto lutetiumElement = nistManager->FindOrBuildElement("Lu");

    const auto csI = new G4Material("CsI", 4.51 * g/cm3, 3, kStateSolid);
    csI->AddElement(cesiumElement, 0.507556);
    csI->AddElement(iodideElement, 0.484639);
    csI->AddElement(thaliumElement, 0.007805);

    const auto labr = new G4Material("LaBr3", 5.08 * g/cm3, 3, kStateSolid);
    labr->AddElement(bromideElement, 0.631308);
    labr->AddElement(cesiumElement, 0.036582);
    labr->AddElement(lanthanumElement, 0.332110);

    const auto lyso = new G4Material("LYSO", 7.1 * g/cm3, 5, kStateSolid);
    lyso->AddElement(oxygenElement, 0.175801);
    lyso->AddElement(siliconElement, 0.061720);
    lyso->AddElement(yttriumElement, 0.019538);
    lyso->AddElement(lutetiumElement, 0.730562);
    lyso->AddElement(cesiumElement, 0.012379);



    // Define Optical Properties
    std::vector<G4double> energy;
    std::vector<G4double> rindex;
    std::vector<G4double> efficiency;
    std::vector<G4double> reflectivity;
    std::vector<G4double> absoption;

//============================================ silicone oil ============================================

    // silicone oil refractive index
    energy = {8.01532E-07*eV, 1.89386E-06*eV, 1.92915E-06*eV, 2.1093E-06*eV, 2.27541E-06*eV, 2.55633E-06*eV, 2.58828E-06*eV};
    rindex = {1.3912, 1.3992, 1.3997, 1.4015, 1.4034, 1.4071, 1.4076};
    const auto siliconeOilRefIndex = new G4MaterialPropertyVector(&energy.front(), &rindex.front(), energy.size());

    // silicone oil absoption length
    energy = {1. * eV, 20. * eV};
    absoption = {15. * cm, 15. * cm};
    const auto siliconeOilAbsLength = new G4MaterialPropertyVector(&energy.front(), &absoption.front(), energy.size());

    // add
    auto siliconeOilPropertiesTable = new G4MaterialPropertiesTable();
    siliconeOilPropertiesTable->AddProperty("RINDEX", siliconeOilRefIndex);
    siliconeOilPropertiesTable->AddProperty("ABSLENGTH", siliconeOilAbsLength);
    siliconeOil->SetMaterialPropertiesTable(siliconeOilPropertiesTable);

//============================================ SiPM ============================================

    // // SiPM refractive index
    // energy = {1.587 * eV, 3.095 * eV};
    // rindex = {1.403, 1.406};
    // const auto siPMRefIndex = new G4MaterialPropertyVector(&energy.front(), &rindex.front(), energy.size());

    // // SiPM absoption length
    // energy = {1.0 * eV, 20.0 * eV};
    // absoption = {0*cm, 0*cm};
    // const auto siPMAbsLength = new G4MaterialPropertyVector(&energy.front(), &absoption.front(), energy.size());

    // // add
    // auto siPMPropertiesTable = new G4MaterialPropertiesTable();
    // siPMPropertiesTable->AddProperty("RINDEX", siPMRefIndex);
    // siPMPropertiesTable->AddProperty("ABSLENGTH", siPMAbsLength);
    // silicon->SetMaterialPropertiesTable(siPMPropertiesTable);

//============================================ Crystal ============================================

    //CsI(Tl)
    auto csiPropertiesTable = new G4MaterialPropertiesTable;
    auto csiProperties(CreateMapFromCSV<G4double>("data/CsI_properties.csv"));
    csiPropertiesTable->AddProperty(
        "RINDEX",
        &csiProperties["energy"].front(),
        &csiProperties["RINDEX"].front(),
        csiProperties["RINDEX"].size()
    );
    csiPropertiesTable->AddProperty(
    "ABSLENGTH",
    &csiProperties["energy"].front(),
    &csiProperties["ABSLENGTH"].front(),
    csiProperties["ABSLENGTH"].size()
    );
    csiPropertiesTable->AddProperty(
        "SCINTILLATIONCOMPONENT1",
        &csiProperties["energy"].front(),
        &csiProperties["SCINTILLATIONCOMPONENT1"].front(),
        csiProperties["SCINTILLATIONCOMPONENT1"].size()
    );
    csiPropertiesTable->AddConstProperty("SCINTILLATIONYIELD", csiProperties["SCINTILLATIONYIELD"].front());
    csiPropertiesTable->AddConstProperty("SCINTILLATIONTIMECONSTANT1", csiProperties["SCINTILLATIONTIMECONSTANT1"].front());
    csiPropertiesTable->AddConstProperty("RESOLUTIONSCALE", csiProperties["RESOLUTIONSCALE"].front());
    csI->SetMaterialPropertiesTable(csiPropertiesTable);

    //LaBr3(Ce)
    auto labrPropertiesTable = new G4MaterialPropertiesTable;
    auto labrProperties(CreateMapFromCSV<G4double>("data/LaBr_properties.csv"));
    labrPropertiesTable->AddProperty(
        "RINDEX",
        &labrProperties["energy"].front(),
        &labrProperties["RINDEX"].front(),
        labrProperties["RINDEX"].size()
    );
    labrPropertiesTable->AddProperty(
        "ABSLENGTH",
        &labrProperties["energy"].front(),
        &labrProperties["ABSLENGTH"].front(),
        labrProperties["ABSLENGTH"].size()
    );
    labrPropertiesTable->AddProperty(
        "SCINTILLATIONCOMPONENT1",
        &labrProperties["energy"].front(),
        &labrProperties["SCINTILLATIONCOMPONENT1"].front(),
        labrProperties["SCINTILLATIONCOMPONENT1"].size()
    );
    labrPropertiesTable->AddConstProperty("SCINTILLATIONYIELD", labrProperties["SCINTILLATIONYIELD"].front());
    labrPropertiesTable->AddConstProperty("SCINTILLATIONTIMECONSTANT1", labrProperties["SCINTILLATIONTIMECONSTANT1"].front());
    labrPropertiesTable->AddConstProperty("RESOLUTIONSCALE", labrProperties["RESOLUTIONSCALE"].front());
    labr->SetMaterialPropertiesTable(labrPropertiesTable);

    //LYSO(Ce)
    auto lysoPropertiesTable = new G4MaterialPropertiesTable;
    auto lysoProperties(CreateMapFromCSV<G4double>("data/LYSO_properties.csv"));
    lysoPropertiesTable->AddProperty(
        "RINDEX",
        &lysoProperties["energy"].front(),
        &lysoProperties["RINDEX"].front(),
        lysoProperties["RINDEX"].size()
    );
    lysoPropertiesTable->AddProperty(
        "ABSLENGTH",
        &lysoProperties["energy"].front(),
        &lysoProperties["ABSLENGTH"].front(),
        lysoProperties["ABSLENGTH"].size()
    );
    lysoPropertiesTable->AddProperty(
        "SCINTILLATIONCOMPONENT1",
        &lysoProperties["energy"].front(),
        &lysoProperties["SCINTILLATIONCOMPONENT1"].front(),
        lysoProperties["SCINTILLATIONCOMPONENT1"].size()
    );
    lysoPropertiesTable->AddConstProperty("SCINTILLATIONYIELD", lysoProperties["SCINTILLATIONYIELD"].front());
    lysoPropertiesTable->AddConstProperty("SCINTILLATIONTIMECONSTANT1", lysoProperties["SCINTILLATIONTIMECONSTANT1"].front());
    lysoPropertiesTable->AddConstProperty("RESOLUTIONSCALE", lysoProperties["RESOLUTIONSCALE"].front());
    lyso->SetMaterialPropertiesTable(lysoPropertiesTable);



    // Define World Volume

    const auto worldSize = 2 * m;
    const auto worldS = new G4Box("world", 0.5 * worldSize, 0.5 * worldSize, 0.5 * worldSize);
    const auto worldLV = new G4LogicalVolume(worldS, galactic, "World");
    const auto worldPV = new G4PVPlacement(
        nullptr,
        {},
        worldLV,
        "World",
        nullptr,
        false,
        0,
        true);

    const auto fCrystalHalfWidth = 7.5*cm;
    const auto fCrystalCoatSpacing = 10.*um;
    const auto fReflectiveCoatThickenss = 20. * um;
    const auto fProtectiveCoatThickenss = 50. * um;
    const auto fSiPMThickness = 250.*um;
    const auto translation = G4ThreeVector(0., 0., fCrystalHalfWidth+fSiPMThickness);

    const auto albox = new G4Box("albox", fCrystalHalfWidth+fReflectiveCoatThickenss, fCrystalHalfWidth+fReflectiveCoatThickenss, fCrystalHalfWidth+fReflectiveCoatThickenss);

    const auto crystalSV = new G4Box("crystal", fCrystalHalfWidth, fCrystalHalfWidth, fCrystalHalfWidth);
    const auto crystalLV = new G4LogicalVolume(crystalSV, lyso, "crystal");
    new G4PVPlacement(
    nullptr,
    G4ThreeVector(),
    crystalLV,
    "crystal",
    worldLV,
    false,
    0,
    true);

    const auto sipmSV = new G4Box("sipm", 3.*mm, 3.*mm, 250.*um);
    const auto sipmLV = new G4LogicalVolume(sipmSV, silicon, "sipm");
    new G4PVPlacement(
    nullptr,
    translation,
    sipmLV,
    "sipm",
    worldLV,
    false,
    0,
    true);


    const auto unionSolid = new G4UnionSolid("unionSolid", crystalSV, sipmSV, 0, translation);

    const auto alSV = new G4SubtractionSolid("al", albox, unionSolid);
    const auto alLV = new G4LogicalVolume(alSV, aluminum, "al");
    new G4PVPlacement(
    nullptr,
    G4ThreeVector(),
    alLV,
    "al",
    worldLV,
    false,
    0,
    true);


    // Define Surface

    energy = {1.0*eV, 20.0*eV};
    reflectivity = {0.95, 0.95};

    auto alSurfacePropertiesTable = new G4MaterialPropertiesTable();
    auto alSurface = new G4OpticalSurface("Al foil", unified, polished, dielectric_metal);
    new G4LogicalSkinSurface("AlSkinSurface", logicReflectiveCoat, alSurface);
    alSurfacePropertiesTable->AddProperty("REFLECTIVITY", energy, reflectivity);
    alSurface->SetMaterialPropertiesTable(alSurfacePropertiesTable);

    energy = {1.0*eV, 20.0*eV};
    reflectivity = {0, 0};
    efficiency = {1, 1}

    auto sipmSurfacePropertiesTable = new G4MaterialPropertiesTable();
    auto sipmSurface = new G4OpticalSurface("SiPM", unified, polished, dielectric_metal);
    new G4LogicalSkinSurface("sipmSkinSurface", logicSiPM, sipmSurface);
    sipmSurfacePropertiesTable->AddProperty("REFLECTIVITY", energy, reflectivity);
    sipmSurfacePropertiesTable->AddProperty("EFFICIENCY", energy, reflectivity);
    sipmSurface->SetMaterialPropertiesTable(sipmSurfacePropertiesTable);


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