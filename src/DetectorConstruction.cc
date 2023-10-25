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

    // O 8
    const auto o16 = new G4Isotope("O16", 8, 16, 15.99*g/mole);
    const auto o17 = new G4Isotope("O17", 8, 17, 16.999*g/mole);
    const auto o18 = new G4Isotope("O18", 8, 18, 17.999*g/mole);
    const auto elo16 = new G4Element("O16", "016", 1);
    elo16->AddIsotope(o16, 1);
    const auto elo17 = new G4Element("O17", "017", 1);
    elo17->AddIsotope(o17, 1);
    const auto elo18 = new G4Element("O18", "018", 1);
    elo18->AddIsotope(o18, 1);

    // Br 35
    const auto br79 = new G4Isotope("Br79", 35, 79, 78.9*g/mole);
    const auto br81 = new G4Isotope("Br81", 35, 81, 80.9*g/mole);
    const auto elbr79 = new G4Element("Br79", "Br79", 1);
    elbr79->AddIsotope(br79, 1);
    const auto elbr81 = new G4Element("Br81", "Br81", 1);
    elbr81->AddIsotope(br81, 1);

    // Y 39
    const auto y89 = new G4Isotope("Y89", 39, 89, 88.9*g/mole);
    const auto ely89 = new G4Element("Y89", "Y89", 1);
    ely89->AddIsotope(y89, 1);

    // I 53
    const auto i127 = new G4Isotope("I127", 53, 127, 126.9*g/mole);
    const auto eli127 = new G4Element("I127", "I127", 1);
    eli127->AddIsotope(i127, 1);

    // Cs 55
    const auto cs133 = new G4Isotope("Cs133", 55, 133, 132.9*g/mole);
    const auto elcs133 = new G4Element("Cs133", "Cs133", 1);
    elcs133->AddIsotope(cs133, 1);

    // La 57
    const auto la138 = new G4Isotope("La138", 57, 138, 137.9*g/mole);
    const auto la139 = new G4Isotope("La139", 57, 139, 138.9*g/mole);
    const auto ella138 = new G4Element("La138", "La138", 1);
    ella138->AddIsotope(la138, 1);
    const auto ella139 = new G4Element("La139", "La139", 1);
    ella139->AddIsotope(la139, 1);

    // Ce 58
    const auto ce136 = new G4Isotope("Ce136", 58, 136, 135.9*g/mole);
    const auto ce138 = new G4Isotope("Ce138", 58, 138, 137.9*g/mole);
    const auto ce140 = new G4Isotope("Ce140", 58, 140, 139.9*g/mole);
    const auto ce142 = new G4Isotope("Ce142", 58, 142, 141.9*g/mole);
    const auto elce136 = new G4Element("Ce136", "Ce136", 1);
    elce136->AddIsotope(ce136, 1);
    const auto elce138 = new G4Element("Ce138", "Ce138", 1);
    elce138->AddIsotope(ce138, 1);
    const auto elce140 = new G4Element("Ce140", "Ce140", 1);
    elce140->AddIsotope(ce140, 1);
    const auto elce142 = new G4Element("Ce142", "Ce142", 1);
    elce142->AddIsotope(ce142, 1);

    // Lu 71
    const auto lu175 = new G4Isotope("Lu175", 71, 175, 174.9*g/mole);
    const auto lu176 = new G4Isotope("Lu176", 71, 176, 175.9*g/mole);
    const auto ellu175 = new G4Element("Lu175", "Lu175", 1);
    ellu175->AddIsotope(lu175, 1);
    const auto ellu176 = new G4Element("Lu176", "Lu176", 1);
    ellu176->AddIsotope(lu176, 1);

    // Tl 81
    const auto tl203 = new G4Isotope("Tl203", 81, 203, 202.97*g/mole);
    const auto tl205 = new G4Isotope("Tl205", 81, 205, 204.97*g/mole);
    const auto eltl203 = new G4Element("Tl203", "Tl203", 1);
    eltl203->AddIsotope(tl203, 1);
    const auto eltl205 = new G4Element("Tl205", "Tl205", 1);
    eltl205->AddIsotope(tl205, 1);

    // const auto csI = new G4Material("CsI", 4.51 * g/cm3, 4, kStateSolid);
    // csI->AddElement(elcs133, 0.507556);
    // csI->AddElement(eli127, 0.484639);
    // csI->AddElement(eltl203, 0.002288);
    // csI->AddElement(eltl205, 0.005517);

    // const auto labr = new G4Material("LaBr3", 5.08 * g/cm3, 8, kStateSolid);
    // labr->AddElement(elbr79, 0.316012);
    // labr->AddElement(elbr81, 0.315192);
    // labr->AddElement(elce136, 0.000033);
    // labr->AddElement(elce138, 0.000045);
    // labr->AddElement(elce140, 0.016151);
    // labr->AddElement(elce142, 0.002059);
    // labr->AddElement(ella138, 0.000309);
    // labr->AddElement(ella139, 0.350195);

    // const auto lyso = new G4Material("LYSO", 7.1 * g/cm3, 11, kStateSolid);
    // lyso->AddElement(elo16, 0.175325);
    // lyso->AddElement(elo17, 0.000071);
    // lyso->AddElement(elo18, 0.000405);
    // lyso->AddElement(siliconElement, 0.061720);
    // lyso->AddElement(ely89, 0.019538);
    // lyso->AddElement(ellu175, 0.711469);
    // lyso->AddElement(ellu176, 0.019093);
    // lyso->AddElement(elce136, 0.000022);
    // lyso->AddElement(elce138, 0.000031);
    // lyso->AddElement(elce140, 0.010932);
    // lyso->AddElement(elce142, 0.001393);

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
    labr->AddElement(bromideElement, 0.631208);
    labr->AddElement(lanthanumElement, 0.350504);
    labr->AddElement(cesiumElement, 0.018288);

    const auto lyso = new G4Material("LYSO", 7.1 * g/cm3, 5, kStateSolid);
    lyso->AddElement(oxygenElement, 0.175801);
    lyso->AddElement(siliconElement, 0.061720);
    lyso->AddElement(yttriumElement, 0.019538);
    lyso->AddElement(lutetiumElement, 0.730562);
    lyso->AddElement(cesiumElement, 0.012379);



    // Define Optical Properties
    std::vector<G4double> energy;
    std::vector<G4double> rindex;
    std::vector<G4double> reflectivity;
    std::vector<G4double> absoption;

//============================================ air ============================================

    energy = {1.0 * eV, 20.0 * eV};
    rindex = {1.0, 1.0};
    const auto airRefIndex = new G4MaterialPropertyVector(&energy.front(), &rindex.front(), energy.size());
    auto airPropertiesTable = new G4MaterialPropertiesTable();
    airPropertiesTable->AddProperty("RINDEX", airRefIndex);
    galactic->SetMaterialPropertiesTable(airPropertiesTable);

//============================================ Al foil ============================================

    energy = {1.0*eV, 20.0*eV};
    reflectivity = {0.95, 0.95};
    const auto alRefIndex = new G4MaterialPropertyVector(&energy.front(), &reflectivity.front(), energy.size());
    auto alPropertiesTable = new G4MaterialPropertiesTable();
    alPropertiesTable->AddProperty("REFLECTIVITY", alRefIndex);
    aluminum->SetMaterialPropertiesTable(alPropertiesTable);

//============================================ silicone oil ============================================

    // silicone oil refractive index
    energy = {1.587 * eV, 3.095 * eV};
    rindex = {1.403, 1.406};
    const auto siliconeOilRefIndex = new G4MaterialPropertyVector(&energy.front(), &rindex.front(), energy.size());

    // silicone oil absoption length
    energy = {1.587 * eV, 3.095 * eV};
    absoption = {50. * cm, 50. * cm};
    const auto siliconeOilAbsLength = new G4MaterialPropertyVector(&energy.front(), &absoption.front(), energy.size());

    // add
    auto siliconeOilPropertiesTable = new G4MaterialPropertiesTable();
    siliconeOilPropertiesTable->AddProperty("RINDEX", siliconeOilRefIndex);
    siliconeOilPropertiesTable->AddProperty("ABSLENGTH", siliconeOilAbsLength);
    siliconeOil->SetMaterialPropertiesTable(siliconeOilPropertiesTable);

//============================================ SiPM ============================================

    // SiPM refractive index
    energy = {1.587 * eV, 3.095 * eV};
    rindex = {1.403, 1.406};
    const auto siPMRefIndex = new G4MaterialPropertyVector(&energy.front(), &rindex.front(), energy.size());

    // SiPM absoption length
    energy = {1.0 * eV, 20.0 * eV};
    absoption = {0*cm, 0*cm};
    const auto siPMAbsLength = new G4MaterialPropertyVector(&energy.front(), &absoption.front(), energy.size());

    // add
    auto siPMPropertiesTable = new G4MaterialPropertiesTable();
    siPMPropertiesTable->AddProperty("RINDEX", siPMRefIndex);
    siPMPropertiesTable->AddProperty("ABSLENGTH", siPMAbsLength);
    silicon->SetMaterialPropertiesTable(siPMPropertiesTable);

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

    auto alSurface = new G4OpticalSurface("Al foil", unified, polished, dielectric_metal);
    alSurface->SetMaterialPropertiesTable(alPropertiesTable);
    new G4LogicalSkinSurface("AlSkinSurface", alLV, alSurface);

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