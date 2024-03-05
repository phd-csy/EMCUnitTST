#include "DetectorConstruction.hh"

#include "CreateMapFromCSV.hh"
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
#include "G4PVDivision.hh"
#include "G4PVPlacement.hh"
#include "G4QuadrangularFacet.hh"
#include "G4ReplicatedSlice.hh"
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

    const auto air = nistManager->BuildMaterialWithNewDensity("Vacuum", "G4_AIR", 1e-12 * g / cm3);
    const auto aluminum = nistManager->FindOrBuildMaterial("G4_Al");
    const auto silicon = nistManager->FindOrBuildMaterial("G4_Si");
    const auto pvc = nistManager->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
    const auto bgo = nistManager->FindOrBuildMaterial("G4_BGO");
    const auto glass = nistManager->FindOrBuildMaterial("G4_GLASS_PLATE");

    const auto carbonElement = nistManager->FindOrBuildElement("C");
    const auto hydrogenElement = nistManager->FindOrBuildElement("H");
    const auto oxygenElement = nistManager->FindOrBuildElement("O");
    const auto nitrogenElement = nistManager->FindOrBuildElement("N");
    const auto siliconElement = nistManager->FindOrBuildElement("Si");
    const auto cesiumElement = nistManager->FindOrBuildElement("Cs");
    const auto iodideElement = nistManager->FindOrBuildElement("I");
    const auto thaliumElement = nistManager->FindOrBuildElement("Tl");
    const auto bromideElement = nistManager->FindOrBuildElement("Br");
    const auto lanthanumElement = nistManager->FindOrBuildElement("La");
    const auto yttriumElement = nistManager->FindOrBuildElement("Y");
    const auto lutetiumElement = nistManager->FindOrBuildElement("Lu");

    const auto potassiumElement = nistManager->FindOrBuildElement("K");
    const auto antimonyElement = nistManager->FindOrBuildElement("Sb");
    const auto bismuthElement = nistManager->FindOrBuildElement("Bi");

    const auto pet = new G4Material("PET", 1.38 * g / cm3, 3, kStateSolid);
    pet->AddElement(hydrogenElement, 0.041962);
    pet->AddElement(carbonElement, 0.625008);
    pet->AddElement(oxygenElement, 0.333030);

    const auto siliconeOil = new G4Material("silicone_oil", 0.97 * g / cm3, 4, kStateLiquid);
    siliconeOil->AddElement(carbonElement, 2);
    siliconeOil->AddElement(hydrogenElement, 6);
    siliconeOil->AddElement(oxygenElement, 1);
    siliconeOil->AddElement(siliconElement, 1);

    const auto csI = new G4Material("CsI", 4.51 * g / cm3, 3, kStateSolid);
    csI->AddElement(cesiumElement, 0.507556);
    csI->AddElement(iodideElement, 0.484639);
    csI->AddElement(thaliumElement, 0.007805);

    const auto bialkali = new G4Material("Bialkali", 2.0 * g / cm3, 3, kStateSolid);
    bialkali->AddElement(potassiumElement, 2);
    bialkali->AddElement(cesiumElement, 1);
    bialkali->AddElement(antimonyElement, 1);

    // const auto plastic = new G4Material("EJ200", 1.023 * g / cm3, 2, kStateSolid);
    // plastic->AddElement(hydrogenElement, 0.084744);
    // plastic->AddElement(carbonElement, 0.915256);

    const auto plastic = new G4Material("Doped", 1.03 * g / cm3, 3, kStateSolid);
    plastic->AddElement(carbonElement, 0.915 * 0.8);
    plastic->AddElement(hydrogenElement, 0.085 * 0.8);
    plastic->AddElement(bismuthElement, 0.2);

    const auto epoxy = new G4Material("Epoxy", 1.18 * g / cm3, 3, kStateSolid);
    epoxy->AddElement(carbonElement, 73.62 * perCent);
    epoxy->AddElement(hydrogenElement, 6.75 * perCent);
    epoxy->AddElement(oxygenElement, 19.63 * perCent);

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

    //////////////////////////////////////////////////
    // Construct Material Optical Properties Tables
    //////////////////////////////////////////////////

    constexpr auto fLambda_min = 200 * nm;
    constexpr auto fLambda_max = 700 * nm;
    std::vector<G4double> fEnergyPair = {h_Planck * c_light / fLambda_max,
                                         h_Planck * c_light / fLambda_min};

    //============================================ Air ================================================

    const auto airPropertiesTable = new G4MaterialPropertiesTable();
    airPropertiesTable->AddProperty("RINDEX", fEnergyPair, {1.00, 1.00});
    air->SetMaterialPropertiesTable(airPropertiesTable);

    //============================================ Optical Coupler ====================================

    // std::vector<G4double> couplerEnergyBin = {8.01532E-07, 1.89386E-06, 1.92915E-06, 2.1093E-06, 2.27541E-06, 2.55633E-06, 2.58828E-06};
    // std::vector<G4double> couplerRefractiveIndex = {1.3912, 1.3992, 1.3997, 1.4015, 1.4034, 1.4071, 1.4076};

    const auto siliconeOilPropertiesTable = new G4MaterialPropertiesTable();
    siliconeOilPropertiesTable->AddProperty("RINDEX", fEnergyPair, {1.465, 1.465});
    siliconeOilPropertiesTable->AddProperty("ABSLENGTH", fEnergyPair, {40 * cm, 40 * cm});
    siliconeOil->SetMaterialPropertiesTable(siliconeOilPropertiesTable);

    //============================================ Optical Window =====================================

    const auto windowPropertiesTable = new G4MaterialPropertiesTable();
    windowPropertiesTable->AddProperty("RINDEX", fEnergyPair, {1.55, 1.55});
    epoxy->SetMaterialPropertiesTable(windowPropertiesTable);

    //============================================ Crystal ============================================

    // CsI(Tl)

    const auto csiPropertiesTable = new G4MaterialPropertiesTable();
    auto csiProperties(CreateMapFromCSV<G4double>("data/CsI_properties.csv"));
    csiPropertiesTable->AddProperty("RINDEX", fEnergyPair, {1.79, 1.79});
    csiPropertiesTable->AddProperty("ABSLENGTH", fEnergyPair, {370 * mm, 370 * mm});
    csiPropertiesTable->AddProperty("SCINTILLATIONCOMPONENT1",
                                    csiProperties["energy"],
                                    csiProperties["SCINTILLATIONCOMPONENT1"]);
    csiPropertiesTable->AddConstProperty("SCINTILLATIONYIELD",
                                         csiProperties["SCINTILLATIONYIELD"].front());
    csiPropertiesTable->AddConstProperty("SCINTILLATIONTIMECONSTANT1",
                                         csiProperties["SCINTILLATIONTIMECONSTANT1"].front());
    csiPropertiesTable->AddConstProperty("RESOLUTIONSCALE",
                                         csiProperties["RESOLUTIONSCALE"].front());
    csI->SetMaterialPropertiesTable(csiPropertiesTable);

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

    // BGO

    std::vector<G4double> bgoWaveLengthBin = {635.394, 629.507, 621.889, 614.27, 606.652, 599.033, 591.414, 583.796, 576.177, 568.905, 562.325, 556.092, 550.205, 544.664, 539.124, 533.929, 529.081, 523.886, 517.999, 511.073, 503.455, 495.836, 488.218, 480.599, 472.981, 465.362, 458.436, 452.895, 448.394, 444.584, 440.775, 436.619, 432.81, 429.347, 425.884, 422.768, 419.651, 416.188, 412.725, 409.262, 405.799, 401.99, 397.488, 392.64, 387.445, 381.558, 376.284};
    std::vector<G4double> bgoEnergyBin(bgoWaveLengthBin.size());
    std::vector<G4double> bgoScintillationComponent1 = {0.107, 0.127, 0.156, 0.186, 0.218, 0.253, 0.289, 0.325, 0.363, 0.402, 0.442, 0.481, 0.521, 0.56, 0.6, 0.641, 0.68, 0.721, 0.758, 0.798, 0.839, 0.875, 0.899, 0.909, 0.9, 0.875, 0.839, 0.8, 0.76, 0.721, 0.68, 0.635, 0.592, 0.552, 0.509, 0.469, 0.428, 0.383, 0.339, 0.296, 0.257, 0.217, 0.176, 0.136, 0.099, 0.059, 0.026};
    std::transform(bgoWaveLengthBin.begin(), bgoWaveLengthBin.end(), bgoEnergyBin.begin(),
                   [](auto val) { return h_Planck * c_light / (val * nm / mm); });

    const auto bgoPropertiesTable = new G4MaterialPropertiesTable();
    bgoPropertiesTable->AddProperty("RINDEX", fEnergyPair, {2.15, 2.15});
    bgoPropertiesTable->AddProperty("ABSLENGTH", fEnergyPair, {40 * cm, 40 * cm});
    bgoPropertiesTable->AddProperty("SCINTILLATIONCOMPONENT1", bgoEnergyBin, bgoScintillationComponent1);
    bgoPropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 10000. / MeV);
    bgoPropertiesTable->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 300 * ns);
    bgoPropertiesTable->AddConstProperty("RESOLUTIONSCALE", 1.0);
    bgo->SetMaterialPropertiesTable(bgoPropertiesTable);

    // Plastic Scintillators

    std::vector<G4double> ej200WaveLengthBin = {499.837, 497.667, 495.162, 492.656, 490.15,
                                                487.645, 485.139, 482.633, 480.128, 477.622,
                                                475.116, 472.611, 470.105, 467.6, 465.094,
                                                462.588, 460.083, 457.577, 455.071, 452.566,
                                                450.06, 447.554, 445.205, 443.013, 440.977,
                                                439.254, 437.532, 435.652, 433.617, 431.268,
                                                429.2, 426.569, 424.064, 422.644, 421.871,
                                                420.587, 419.71, 419.052, 418.395, 417.518,
                                                416.86, 416.202, 415.294, 414.386, 413.697,
                                                412.788, 411.849, 410.752, 409.813, 408.717,
                                                407.464, 406.211, 404.802, 402.979, 400.26,
                                                397.754};
    std::vector<G4double> ej200EnergyBin(ej200WaveLengthBin.size());
    std::vector<G4double> ej200ScintillationComponent1 = {0.06, 0.063, 0.074, 0.086, 0.099,
                                                          0.112, 0.13, 0.15, 0.173, 0.198,
                                                          0.228, 0.263, 0.306, 0.348, 0.382,
                                                          0.41, 0.431, 0.455, 0.483, 0.514,
                                                          0.55, 0.592, 0.637, 0.685, 0.731,
                                                          0.773, 0.816, 0.86, 0.9, 0.939,
                                                          0.979, 0.998, 0.997, 0.979, 0.951,
                                                          0.898, 0.846, 0.8, 0.754, 0.698,
                                                          0.651, 0.602, 0.545, 0.488, 0.438,
                                                          0.385, 0.339, 0.292, 0.25, 0.207,
                                                          0.161, 0.117, 0.072, 0.026, 0.011,
                                                          0.003};

    std::transform(ej200WaveLengthBin.begin(), ej200WaveLengthBin.end(), ej200EnergyBin.begin(),
                   [](auto val) { return h_Planck * c_light / (val * nm / mm); });

    const auto plasticPropertiesTable = new G4MaterialPropertiesTable();
    plasticPropertiesTable->AddProperty("RINDEX", fEnergyPair, {1.58, 1.58});
    plasticPropertiesTable->AddProperty("ABSLENGTH", fEnergyPair, {40 * cm, 40 * cm});
    plasticPropertiesTable->AddProperty("SCINTILLATIONCOMPONENT1", ej200EnergyBin, ej200ScintillationComponent1);
    plasticPropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 20000. / MeV);
    plasticPropertiesTable->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 3 * ns);
    plasticPropertiesTable->AddConstProperty("RESOLUTIONSCALE", 1.0);
    plastic->SetMaterialPropertiesTable(plasticPropertiesTable);

    //============================================ Surface ============================================

    const auto rfSurfacePropertiesTable = new G4MaterialPropertiesTable();
    rfSurfacePropertiesTable->AddProperty("REFLECTIVITY", fEnergyPair, {0.985, 0.985});

    const auto couplerSurfacePropertiesTable = new G4MaterialPropertiesTable();
    couplerSurfacePropertiesTable->AddProperty("TRANSMITTANCE", fEnergyPair, {1, 1});

    const auto airPaintSurfacePropertiesTable = new G4MaterialPropertiesTable();
    airPaintSurfacePropertiesTable->AddProperty("REFLECTIVITY", fEnergyPair, {0., 0.});

    const auto cathodeSurfacePropertiesTable = new G4MaterialPropertiesTable();

    std::vector<G4double> sipmEnergyBin = {1.40E-06, 1.44E-06, 1.50E-06, 1.56E-06, 1.63E-06, 1.69E-06, 1.76E-06, 1.82E-06,
                                           1.88E-06, 1.93E-06, 1.98E-06, 2.03E-06, 2.09E-06, 2.17E-06, 2.24E-06, 2.30E-06,
                                           2.38E-06, 2.48E-06, 2.64E-06, 2.84E-06, 3.01E-06, 3.13E-06, 3.22E-06, 3.25E-06,
                                           3.34E-06, 3.41E-06, 3.49E-06, 3.59E-06, 3.62E-06, 3.67E-06, 3.72E-06, 3.77E-06,
                                           3.79E-06, 3.86E-06};
    std::vector<G4double> sipmQuantumEfficiency = {0.038204608, 0.049561724, 0.066092061, 0.084585974, 0.103160774, 0.12483983, 0.145920339, 0.167429879,
                                                   0.187343302, 0.206777439, 0.229192371, 0.249945581, 0.269091999, 0.297374916, 0.323900966, 0.344069886,
                                                   0.364481152, 0.385104958, 0.399383156, 0.392068233, 0.369892386, 0.348694706, 0.324319675, 0.301670517,
                                                   0.279767001, 0.258464601, 0.231195542, 0.198458795, 0.162145893, 0.13310099, 0.098584057, 0.072563977,
                                                   0.053495328, 0.032820352};

    cathodeSurfacePropertiesTable->AddProperty("REFLECTIVITY", fEnergyPair, {0., 0.});
    cathodeSurfacePropertiesTable->AddProperty("EFFICIENCY", sipmEnergyBin, sipmQuantumEfficiency);

    /////////////////////////////////////////////
    // Construct Volumes
    /////////////////////////////////////////////

    const auto worldSize = 2 * m;
    const auto worldS = new G4Box("world", 0.5 * worldSize, 0.5 * worldSize, 0.5 * worldSize);
    const auto worldLV = new G4LogicalVolume(worldS, air, "World");
    const auto worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World", nullptr, false, 0, true);

    const auto fCrystalWidth = 3 * cm;
    const auto fCrystalLength = 8 * cm;

    const auto fReflectorThickness = 0.5 * mm;

    const auto fCouplerThickness = 0.1 * mm;
    const auto fWindowThickness = 1 * mm;
    const auto fSiPMWidth = 3 * mm;
    const auto fSiPMThickness = 0.1 * mm;

    // const auto nsect = 5;
    // const auto fCrystalLength = 146.73 * mm;
    // const auto fPolygonCrystalWidth = 20.522 * mm;

    // std::vector<G4TwoVector> polygon(nsect);
    // G4double ang = twopi / nsect;
    // for (int i = 0; i < nsect; ++i) {
    //     G4double phi = i * ang;
    //     G4double cosphi = std::cos(phi);
    //     G4double sinphi = std::sin(phi);
    //     polygon[i].set(fPolygonCrystalWidth * cosphi, fPolygonCrystalWidth * sinphi);
    // }

    // G4TwoVector offsetA(0, 0), offsetB(0, offsetlength);
    // G4double scaleA = 1, scaleB = 1 + (100 / (100 + fCrystalLength));
    // G4double scaleA = 1, scaleB = 2;

    // const auto fPMTRadius = fPolygonCrystalWidth * scaleB * 0.5;
    // const auto fPMTRadius = 25.5 * mm;
    // const auto fPMTCathodeRadius = 23 * mm;

    const auto Transform =
        [&fCrystalLength, crystalTail = G4ThreeVector(0, 0, fCrystalLength / 2)](double transformDistance) {
            return G4Translate3D{crystalTail + G4ThreeVector(0, 0, transformDistance)};
        };

    // const auto crystalSV = new G4ExtrudedSolid("Extruded", polygon, fCrystalLength / 2, offsetA, scaleA, offsetB, scaleB);

    const auto crystalSV = new G4Box("crystal", fCrystalWidth / 2, fCrystalWidth / 2, fCrystalLength / 2);
    //========================================== CsI(Tl) ============================================
    const auto crystalLV = new G4LogicalVolume{crystalSV, csI, "crystal"};
    //========================================== LaBr3(Ce) ==========================================
    // const auto crystalLV = new G4LogicalVolume{crystalSV, labr, "crystal"};
    //========================================== LYSO(Ce) ===========================================
    // const auto crystalLV = new G4LogicalVolume{crystalSV, lyso, "crystal"};
    //===============================================================================================
    // const auto crystalLV = new G4LogicalVolume{crystalSV, bgo, "crystal"};
    const auto crystalPV = new G4PVPlacement(G4Transform3D{}, crystalLV, "crystal", worldLV, false, 0, true);

    const auto reflectorSV = new G4Box("reflector", fCrystalWidth / 2 + fReflectorThickness, fCrystalWidth / 2 + fReflectorThickness, fCrystalLength / 2 + fReflectorThickness);
    const auto cuttedReflector = new G4SubtractionSolid("reflector", reflectorSV, crystalSV, G4Translate3D(0, 0, fReflectorThickness));
    const auto reflectorLV = new G4LogicalVolume{cuttedReflector, pet, "reflector"};
    const auto reflectorPV = new G4PVPlacement{G4Translate3D(0, 0, -fReflectorThickness), reflectorLV, "reflector", worldLV, false, 0, true};

    const auto solidCoupler = new G4Box("coupler", fCrystalWidth / 2, fCrystalWidth / 2, fCouplerThickness / 2);
    const auto logicalCoupler = new G4LogicalVolume(solidCoupler, siliconeOil, "CrystalCoupler");
    const auto physicalCoupler = new G4PVPlacement(Transform(fCouplerThickness / 2),
                                                   logicalCoupler, "CrystalCoupler", worldLV, false, 0, true);

    const auto windowSV = new G4Box("window", fCrystalWidth / 2, fCrystalWidth / 2, fWindowThickness / 2);
    const auto windowLV = new G4LogicalVolume(windowSV, epoxy, "window");
    const auto windowPV = new G4PVPlacement{Transform(fCouplerThickness + fWindowThickness / 2),
                                            windowLV, "window", worldLV, false, 0, true};

    const auto solidSipm = new G4Box("sipm", fSiPMWidth / 2, fSiPMWidth / 2, fSiPMThickness / 2);
    const auto logicalSipm = new G4LogicalVolume(solidSipm, silicon, "sipm");
    int n = 4;
    new G4PVPlacement(0, G4ThreeVector((n - 0.5) * (fSiPMWidth + 0.2 * mm), (n - 0.5) * (fSiPMWidth + 0.2 * mm), fCrystalLength / 2 + fCouplerThickness + fWindowThickness + fSiPMThickness / 2),
                      logicalSipm, "sipm", worldLV, false, 0, true);

    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            new G4PVPlacement(0, G4ThreeVector((n - 0.5) * (fSiPMWidth + 0.2 * mm) - i * (fSiPMWidth + 0.2 * mm), (n - 0.5) * (fSiPMWidth + 0.2 * mm) - j * (fSiPMWidth + 0.2 * mm), fCrystalLength / 2 + fCouplerThickness + fWindowThickness + fSiPMThickness / 2),
                              logicalSipm, "sipm", worldLV, false, 0, true);
        }
    }

    // Define Surface

    // const auto rfSurface = new G4OpticalSurface("reflector", LUT, polishedvm2000glue, dielectric_LUT);
    const auto rfSurface = new G4OpticalSurface("reflector", unified, polished, dielectric_metal);
    rfSurface->SetMaterialPropertiesTable(rfSurfacePropertiesTable);
    new G4LogicalBorderSurface("rfSurface", crystalPV, reflectorPV, rfSurface);

    const auto couplerSurface = new G4OpticalSurface("coupler", unified, polished, dielectric_dielectric);
    couplerSurface->SetMaterialPropertiesTable(couplerSurfacePropertiesTable);
    new G4LogicalBorderSurface("couplerSurface", crystalPV, physicalCoupler, couplerSurface);

    const auto airPaintSurface = new G4OpticalSurface("Paint", unified, polished, dielectric_metal);
    airPaintSurface->SetMaterialPropertiesTable(airPaintSurfacePropertiesTable);
    new G4LogicalBorderSurface("AirPaintSurface", reflectorPV, crystalPV, airPaintSurface);

    const auto cathodeSurface = new G4OpticalSurface("Cathode", unified, polished, dielectric_metal);
    cathodeSurface->SetMaterialPropertiesTable(cathodeSurfacePropertiesTable);
    new G4LogicalSkinSurface("cathodeSkinSurface", logicalSipm, cathodeSurface);

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
