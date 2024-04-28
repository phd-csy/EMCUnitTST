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

    const auto xenonElement = nistManager->FindOrBuildElement("Xe");

    const auto liquidXe = new G4Material("LXe", 3.02 * g / cm3, 1, kStateLiquid);
    liquidXe->AddElement(xenonElement, 1);

    const auto pet = new G4Material("PET", 1.38 * g / cm3, 3, kStateSolid);
    pet->AddElement(hydrogenElement, 0.041962);
    pet->AddElement(carbonElement, 0.625008);
    pet->AddElement(oxygenElement, 0.333030);

    const auto acrylic = new G4Material("acrylic", 1.19 * g / cm3, 3);
    acrylic->AddElement(carbonElement, 5);
    acrylic->AddElement(hydrogenElement, 8);
    acrylic->AddElement(oxygenElement, 2);

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
    windowPropertiesTable->AddProperty("RINDEX", fEnergyPair, {1.49, 1.49});
    glass->SetMaterialPropertiesTable(windowPropertiesTable);

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

    std::vector<G4double> xeEnergyBin = {7.0 * eV, 7.07 * eV, 7.14 * eV};
    const auto xePropertiesTable = new G4MaterialPropertiesTable();
    xePropertiesTable->AddProperty("RINDEX", xeEnergyBin, {1.59, 1.57, 1.54});
    xePropertiesTable->AddProperty("ABSLENGTH", xeEnergyBin, {35. * cm, 35. * cm, 35. * cm});
    xePropertiesTable->AddProperty("SCINTILLATIONCOMPONENT1", xeEnergyBin, {0.1, 1.0, 0.1});
    xePropertiesTable->AddProperty("SCINTILLATIONCOMPONENT2", xeEnergyBin, {0.1, 1.0, 0.1});
    xePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 45000. / MeV);
    xePropertiesTable->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 20 * ns);
    xePropertiesTable->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 45 * ns);
    xePropertiesTable->AddConstProperty("RESOLUTIONSCALE", 1.0);
    liquidXe->SetMaterialPropertiesTable(xePropertiesTable);

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

    std::vector<G4double> pmtWaveLengthBin = {715.759, 704.541, 687.714, 670.887, 654.06, 637.234, 620.807, 606.384,
                                              592.562, 584.019, 577.939, 571.814, 566.671, 562.542, 558.307, 553.099,
                                              547.49, 541.48, 534.669, 527.057, 519.361, 511.903, 505.422, 499.413,
                                              493.804, 487.821, 481.651, 473.856, 465.246, 456.513, 443.724, 427.297,
                                              410.47, 393.643, 376.816, 359.989, 347.57, 341.159, 335.766, 332.51,
                                              328.633, 325.763, 323.291, 320.44, 318.552, 316.506, 313.615, 312.091,
                                              309.509, 307.334, 305.549, 302.698, 301.596, 299.342, 296.555, 294.161,
                                              291.88, 288.847, 285.737, 281.108};
    std::vector<G4double> pmtQuantumEfficiency = {0.206, 0.237, 0.282, 0.398, 0.74, 1.321, 2.113, 3.024,
                                                  3.969, 4.905, 5.856, 6.819, 7.723, 8.636, 9.57, 10.508,
                                                  11.467, 12.374, 13.281, 14.205, 15.191, 16.195, 17.195,
                                                  18.114, 18.987, 19.886, 20.83, 21.794, 22.8, 23.806, 24.644,
                                                  25.312, 25.713, 25.932, 25.835, 25.279, 24.266, 23.367, 22.357,
                                                  21.43, 20.344, 19.319, 18.363, 17.294, 16.265, 15.232, 14.053,
                                                  12.759, 11.486, 10.345, 9.229, 8.193, 7.198, 6.108, 5.136, 4.241,
                                                  3.37, 2.403, 1.447, 0.466}; // ET 9269B

    // std::vector<double> pmtWaveLengthBin = {724.051, 718.527, 710.242, 701.266, 694.361, 685.386, 677.1, 668.815,
    //                                         658.458, 652.244, 642.578, 632.221, 621.864, 610.127, 596.318, 581.818,
    //                                         567.319, 559.724, 547.986, 538.32, 523.82, 512.083, 495.512, 483.084,
    //                                         461.68, 447.871, 414.039, 398.159, 380.207, 360.184, 338.78, 328.423,
    //                                         317.376, 310.472, 303.567, 298.734};

    // std::vector<double> pmtQuantumEfficiency = {0.011, 0.015, 0.022, 0.033, 0.05, 0.076, 0.115, 0.173,
    //                                             0.289, 0.393, 0.593, 0.899, 1.337, 2.002, 2.956, 4.335,
    //                                             5.978, 6.901, 8.648, 9.915, 12.34, 13.956, 16.558, 18.472,
    //                                             21.179, 23.148, 26.18, 27.841, 27.841, 26.359, 22.678, 19.511,
    //                                             15.15, 12.769, 8.887, 5.286}; // CR 284

    std::vector<G4double> cathodeSurfacePropertiesEnergy(pmtWaveLengthBin.size());
    std::vector<G4double> cathodeSurfacePropertiesEfficiency(pmtQuantumEfficiency.size());

    std::transform(pmtWaveLengthBin.begin(), pmtWaveLengthBin.end(), cathodeSurfacePropertiesEnergy.begin(),
                   [](auto val) { return h_Planck * c_light / (val * nm / mm); });
    std::transform(pmtQuantumEfficiency.begin(), pmtQuantumEfficiency.end(), cathodeSurfacePropertiesEfficiency.begin(),
                   [](auto n) { return n * perCent; });

    cathodeSurfacePropertiesTable->AddProperty("REFLECTIVITY", fEnergyPair, {0., 0.});
    // cathodeSurfacePropertiesTable->AddProperty("EFFICIENCY", cathodeSurfacePropertiesEnergy, cathodeSurfacePropertiesEfficiency);
    cathodeSurfacePropertiesTable->AddProperty("EFFICIENCY", fEnergyPair, {1., 1.});

    /////////////////////////////////////////////
    // Construct Volumes
    /////////////////////////////////////////////

    const auto worldSize = 2 * m;
    const auto worldS = new G4Box("world", 0.5 * worldSize, 0.5 * worldSize, 0.5 * worldSize);
    const auto worldLV = new G4LogicalVolume(worldS, air, "World");
    const auto worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World", nullptr, false, 0, true);

    const auto fCrystalWidth = 3 * cm;
    const auto fCrystalLength = 8 * cm;

    const auto fReflectorThickness = 5 * cm;

    const auto fPMTRadius = 19 * mm;
    const auto fPMTCathodeRadius = 17 * mm;

    const auto fCouplerThickness = 0.1 * mm;
    const auto fWindowThickness = 1 * mm;
    const auto fCathodeThickness = 20 * nm;

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
    // const auto crystalLV = new G4LogicalVolume{crystalSV, csI, "crystal"};
    //========================================== LaBr3(Ce) ==========================================
    // const auto crystalLV = new G4LogicalVolume{crystalSV, labr, "crystal"};
    //========================================== LYSO(Ce) ===========================================
    // const auto crystalLV = new G4LogicalVolume{crystalSV, lyso, "crystal"};
    //===============================================================================================
    const auto crystalLV = new G4LogicalVolume{crystalSV, liquidXe, "crystal"};
    const auto crystalPV = new G4PVPlacement(G4Transform3D{}, crystalLV, "crystal", worldLV, false, 0, true);

    const auto reflectorSV = new G4Box("reflector", fCrystalWidth / 2 + fReflectorThickness, fCrystalWidth / 2 + fReflectorThickness, fCrystalLength / 2 + fReflectorThickness);
    const auto cuttedReflector = new G4SubtractionSolid("reflector", reflectorSV, crystalSV, G4Translate3D(0, 0, fReflectorThickness));
    const auto reflectorLV = new G4LogicalVolume{cuttedReflector, aluminum, "reflector"};
    const auto reflectorPV = new G4PVPlacement{G4Translate3D(0, 0, -fReflectorThickness), reflectorLV, "reflector", worldLV, false, 0, true};

    const auto couplerSV = new G4Tubs("coupler", 0, fPMTRadius, fCouplerThickness / 2, 0, 2 * pi);
    const auto couplerLV = new G4LogicalVolume(couplerSV, siliconeOil, "CrystalCoupler");
    const auto couplerPV = new G4PVPlacement(Transform(fCouplerThickness / 2),
                                             couplerLV, "CrystalCoupler", worldLV, false, 0, true);

    const auto windowSV = new G4Tubs("window", 0, fPMTRadius, fWindowThickness / 2, 0, 2 * pi);
    const auto windowLV = new G4LogicalVolume(windowSV, glass, "window");
    const auto windowPV = new G4PVPlacement{Transform(fCouplerThickness + fWindowThickness / 2),
                                            windowLV, "window", worldLV, false, 0, true};

    const auto sipmSV = new G4Tubs("sipm", 0, fPMTCathodeRadius, fCathodeThickness / 2, 0, 2 * pi);
    const auto sipmLV = new G4LogicalVolume(sipmSV, silicon, "sipm");
    new G4PVPlacement(Transform(fCouplerThickness + fWindowThickness + fCathodeThickness / 2),
                      sipmLV, "sipm", worldLV, false, 0, true);

    // Define Surface

    // const auto rfSurface = new G4OpticalSurface("reflector", LUT, polishedvm2000glue, dielectric_LUT);
    const auto rfSurface = new G4OpticalSurface("reflector", unified, polished, dielectric_metal);
    rfSurface->SetMaterialPropertiesTable(rfSurfacePropertiesTable);
    new G4LogicalBorderSurface("rfSurface", crystalPV, reflectorPV, rfSurface);

    const auto couplerSurface = new G4OpticalSurface("coupler", unified, polished, dielectric_dielectric);
    couplerSurface->SetMaterialPropertiesTable(couplerSurfacePropertiesTable);
    new G4LogicalBorderSurface("couplerSurface", crystalPV, couplerPV, couplerSurface);

    const auto airPaintSurface = new G4OpticalSurface("Paint", unified, polished, dielectric_metal);
    airPaintSurface->SetMaterialPropertiesTable(airPaintSurfacePropertiesTable);
    new G4LogicalBorderSurface("AirPaintSurface", reflectorPV, crystalPV, airPaintSurface);

    const auto cathodeSurface = new G4OpticalSurface("Cathode", unified, polished, dielectric_metal);
    cathodeSurface->SetMaterialPropertiesTable(cathodeSurfacePropertiesTable);
    new G4LogicalSkinSurface("cathodeSkinSurface", sipmLV, cathodeSurface);

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
