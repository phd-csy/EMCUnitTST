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

    const auto potassiumElement = nistManager->FindOrBuildElement("K");
    const auto antimonyElement = nistManager->FindOrBuildElement("Sb");

    const auto siliconeOil = new G4Material("silicone_oil", 0.97 * g / cm3, 4, kStateLiquid);
    siliconeOil->AddElement(carbonElement, 2);
    siliconeOil->AddElement(hydrogenElement, 6);
    siliconeOil->AddElement(oxygenElement, 1);
    siliconeOil->AddElement(siliconElement, 1);

    const auto glass = new G4Material("Fused Silica", 2.64 * g / cm3, 2, kStateSolid);
    glass->AddElement(oxygenElement, 0.532570);
    glass->AddElement(siliconElement, 0.467430);

    const auto acrylic = new G4Material("Acrylic", 1.19 * g / cm3, 3);
    acrylic->AddElement(carbonElement, 5);
    acrylic->AddElement(hydrogenElement, 8);
    acrylic->AddElement(oxygenElement, 2);

    const auto csI = new G4Material("CsI", 4.51 * g / cm3, 3, kStateSolid);
    csI->AddElement(cesiumElement, 0.507556);
    csI->AddElement(iodideElement, 0.484639);
    csI->AddElement(thaliumElement, 0.007805);

    const auto bialkali = new G4Material("Bialkali", 2.0 * g / cm3, 3, kStateSolid);
    bialkali->AddElement(potassiumElement, 2);
    bialkali->AddElement(cesiumElement, 1);
    bialkali->AddElement(antimonyElement, 1);

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

    const auto pmmaPropertiesTable = new G4MaterialPropertiesTable();
    pmmaPropertiesTable->AddProperty("RINDEX", "PMMA");

    constexpr auto fLambda_min = 200 * nm;
    constexpr auto fLambda_max = 700 * nm;
    std::vector<G4double> fEnergyPair = {h_Planck * c_light / fLambda_max,
                                         h_Planck * c_light / fLambda_min};

    //============================================ Air ================================================

    const auto airPropertiesTable = new G4MaterialPropertiesTable();
    airPropertiesTable->AddProperty("RINDEX", fEnergyPair, {1.00, 1.00});
    galactic->SetMaterialPropertiesTable(airPropertiesTable);

    //============================================ Optical Coupler ====================================

    // std::vector<G4double> couplerEnergyBin = {8.01532E-07, 1.89386E-06, 1.92915E-06, 2.1093E-06, 2.27541E-06, 2.55633E-06, 2.58828E-06};
    // std::vector<G4double> couplerRefractiveIndex = {1.3912, 1.3992, 1.3997, 1.4015, 1.4034, 1.4071, 1.4076};

    const auto siliconeOilPropertiesTable = new G4MaterialPropertiesTable();
    siliconeOilPropertiesTable->AddProperty("RINDEX", fEnergyPair, {1.465, 1.465});
    siliconeOilPropertiesTable->AddProperty("ABSLENGTH", fEnergyPair, {15 * cm, 15 * cm});
    siliconeOil->SetMaterialPropertiesTable(siliconeOilPropertiesTable);

    //============================================ Quartz =============================================

    const auto windowPropertiesTable = new G4MaterialPropertiesTable();
    windowPropertiesTable->AddProperty("RINDEX", fEnergyPair, {1.54, 1.54});
    glass->SetMaterialPropertiesTable(windowPropertiesTable);

    //============================================ PMMA ===============================================

    const auto NENTRIES = 11;

    std::vector<G4double> RINDEX_ACRYLIC;
    std::vector<G4double> ENERGY_ACRYLIC;

    // Parameterization for refractive index of High Grade PMMA

    G4double bParam[4] = {1760.7010, -1.3687, 2.4388e-3, -1.5178e-6};
    auto lambda = 0.;

    for (auto i = 0; i < NENTRIES; ++i) {

        // want energies in increasing order
        lambda = fLambda_min + (NENTRIES - 1 - i) * (fLambda_max - fLambda_min) / float(NENTRIES - 1);
        RINDEX_ACRYLIC.push_back(0.0);

        for (auto jj = 0; jj < 4; jj++) {
            RINDEX_ACRYLIC[i] += (bParam[jj] / 1000.0) * std::pow(lambda / nm, jj);
        }

        ENERGY_ACRYLIC.push_back(h_Planck * c_light / lambda); // Convert from wavelength to energy ;
                                                               // G4cout << ENERGY_ACRYLIC[i]/eV << " " << lambda/nm << " " << RINDEX_ACRYLIC[i] << G4endl ;
    }

    const auto acrylicPropertiesTable = new G4MaterialPropertiesTable();
    acrylicPropertiesTable->AddProperty("RINDEX", ENERGY_ACRYLIC, RINDEX_ACRYLIC);

    std::vector<G4double> LAMBDAABS{
        100.0,
        246.528671, 260.605103, 263.853516, 266.019104, 268.726105,
        271.433136, 273.598724, 276.305725, 279.554138, 300.127380,
        320.159241, 340.191101, 360.764343, 381.337585, 399.745239,
        421.401276, 440.891724, 460.382172, 480.414001, 500.987274,
        520.477722, 540.509583, 559.458618,
        700.0};

    std::vector<G4double> ABS // Transmission (in %) of  3mm thick PMMA
        {
            0.0000000,
            0.0000000, 5.295952, 9.657321, 19.937695, 29.283491,
            39.252335, 48.598133, 58.255451, 65.109039, 79.439247,
            85.669785, 89.719627, 91.277260, 91.588783, 91.900307,
            91.588783, 91.277260, 91.277260, 91.588783, 91.588783,
            91.900307, 91.900307, 91.588783,
            91.5};

    acrylicPropertiesTable->AddProperty("ABSLENGTH", new G4MaterialPropertyVector());

    for (size_t i = 0; i < ABS.size(); ++i) {
        auto energy = h_Planck * c_light / (LAMBDAABS[i] * nm);
        auto abslength = 0.0;

        if (ABS[i] <= 0.0) {
            abslength = 1.0 / kInfinity;
        } else {
            abslength = -3.0 * mm / (G4Log(ABS[i] / 100.0));
        }

        acrylicPropertiesTable->AddEntry("ABSLENGTH", energy, abslength);
    }

    acrylic->SetMaterialPropertiesTable(acrylicPropertiesTable);
    acrylicPropertiesTable->DumpTable();

    //============================================ Crystal ============================================

    // CsI(Tl)

    const auto csiPropertiesTable = new G4MaterialPropertiesTable();
    auto csiProperties(CreateMapFromCSV<G4double>("data/CsI_properties.csv"));
    csiPropertiesTable->AddProperty("RINDEX", fEnergyPair, {1.79, 1.79});
    csiPropertiesTable->AddProperty("GROUPVEL", fEnergyPair, {167.482, 167.482});
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

    //============================================ Surface ============================================

    const auto rfSurfacePropertiesTable = new G4MaterialPropertiesTable();
    rfSurfacePropertiesTable->AddProperty("REFLECTIVITY", fEnergyPair, {0.985, 0.985});

    const auto couplerSurfacePropertiesTable = new G4MaterialPropertiesTable();
    couplerSurfacePropertiesTable->AddProperty("TRANSMITTANCE", fEnergyPair, {1, 1});

    const auto airPaintSurfacePropertiesTable = new G4MaterialPropertiesTable();
    airPaintSurfacePropertiesTable->AddProperty("REFLECTIVITY", fEnergyPair, {0., 0.});

    const auto cathodeSurfacePropertiesTable = new G4MaterialPropertiesTable();

    // std::vector<G4double> pmtWaveLengthBin = {715.759, 704.541, 687.714, 670.887, 654.06, 637.234, 620.807, 606.384,
    //                                           592.562, 584.019, 577.939, 571.814, 566.671, 562.542, 558.307, 553.099,
    //                                           547.49, 541.48, 534.669, 527.057, 519.361, 511.903, 505.422, 499.413,
    //                                           493.804, 487.821, 481.651, 473.856, 465.246, 456.513, 443.724, 427.297,
    //                                           410.47, 393.643, 376.816, 359.989, 347.57, 341.159, 335.766, 332.51,
    //                                           328.633, 325.763, 323.291, 320.44, 318.552, 316.506, 313.615, 312.091,
    //                                           309.509, 307.334, 305.549, 302.698, 301.596, 299.342, 296.555, 294.161,
    //                                           291.88, 288.847, 285.737, 281.108};
    // std::vector<G4double> pmtQuantumEfficiency = {0.206, 0.237, 0.282, 0.398, 0.74, 1.321, 2.113, 3.024,
    //                                               3.969, 4.905, 5.856, 6.819, 7.723, 8.636, 9.57, 10.508,
    //                                               11.467, 12.374, 13.281, 14.205, 15.191, 16.195, 17.195,
    //                                               18.114, 18.987, 19.886, 20.83, 21.794, 22.8, 23.806, 24.644,
    //                                               25.312, 25.713, 25.932, 25.835, 25.279, 24.266, 23.367, 22.357,
    //                                               21.43, 20.344, 19.319, 18.363, 17.294, 16.265, 15.232, 14.053,
    //                                               12.759, 11.486, 10.345, 9.229, 8.193, 7.198, 6.108, 5.136, 4.241,
    //                                               3.37, 2.403, 1.447, 0.466};

    // std::vector<G4double> cathodeSurfacePropertiesEnergy(pmtWaveLengthBin.size());
    // std::vector<G4double> cathodeSurfacePropertiesEfficiency(pmtQuantumEfficiency.size());

    // std::transform(pmtWaveLengthBin.begin(), pmtWaveLengthBin.end(), cathodeSurfacePropertiesEnergy.begin(),
    //                [](auto val) { return h_Planck * c_light / (val * nm / mm); });
    // std::transform(pmtQuantumEfficiency.begin(), pmtQuantumEfficiency.end(), cathodeSurfacePropertiesEfficiency.begin(),
    //                [](auto n) { return n * perCent; });

    // cathodeSurfacePropertiesTable->AddProperty("REFLECTIVITY", fEnergyPair, {0., 0.});
    // cathodeSurfacePropertiesTable->AddProperty("EFFICIENCY", cathodeSurfacePropertiesEnergy, cathodeSurfacePropertiesEfficiency);

    // ======================================================================================================================================
    std::vector<double> cathodeEnergyBin = {1.40E-06, 1.44E-06, 1.50E-06, 1.56E-06, 1.63E-06, 1.69E-06, 1.76E-06, 1.82E-06,
                                            1.88E-06, 1.93E-06, 1.98E-06, 2.03E-06, 2.09E-06, 2.17E-06, 2.24E-06, 2.30E-06,
                                            2.38E-06, 2.48E-06, 2.64E-06, 2.84E-06, 3.01E-06, 3.13E-06, 3.22E-06, 3.25E-06,
                                            3.34E-06, 3.41E-06, 3.49E-06, 3.59E-06, 3.62E-06, 3.67E-06, 3.72E-06, 3.77E-06,
                                            3.79E-06, 3.86E-06};

    std::vector<double> cathodeQuantumEfficiency = {0.038204608, 0.049561724, 0.066092061, 0.084585974, 0.103160774, 0.12483983, 0.145920339, 0.167429879,
                                                    0.187343302, 0.206777439, 0.229192371, 0.249945581, 0.269091999, 0.297374916, 0.323900966, 0.344069886,
                                                    0.364481152, 0.385104958, 0.399383156, 0.392068233, 0.369892386, 0.348694706, 0.324319675, 0.301670517,
                                                    0.279767001, 0.258464601, 0.231195542, 0.198458795, 0.162145893, 0.13310099, 0.098584057, 0.072563977,
                                                    0.053495328, 0.032820352};

    cathodeSurfacePropertiesTable->AddProperty("REFLECTIVITY", fEnergyPair, {0., 0.});
    cathodeSurfacePropertiesTable->AddProperty("EFFICIENCY", cathodeEnergyBin, cathodeQuantumEfficiency);

    /////////////////////////////////////////////
    // Construct Volumes
    /////////////////////////////////////////////

    const auto worldSize = 2 * m;
    const auto worldS = new G4Box("world", 0.5 * worldSize, 0.5 * worldSize, 0.5 * worldSize);
    const auto worldLV = new G4LogicalVolume(worldS, galactic, "World");
    const auto worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World", nullptr, false, 0, true);

    const auto fCrystalWidth = 3 * cm;
    const auto fCrystalLength = 8 * cm;

    const auto fHeadCouplerWidth = 3 * cm;
    const auto fTailCouplerWidth = 1 * cm;
    const auto fCouplerThickness = 0.1 * mm;

    const auto fLightGuideHeadWidth = fCrystalWidth;
    const auto fLightGuideTailWidth = 1 * cm;
    const auto fLightGuideThickness = 2.54 * cm;

    const auto fSiPMWidth = 6 * mm;
    const auto fSiPMThickness = 1.45 * mm;

    // const auto offsetlength = 0 * cm;

    // const auto translation1 = G4ThreeVector(0., offsetlength, fCrystalHalfLength + fCouplerHalfThickness);
    // const auto translation2 = G4ThreeVector(0., offsetlength, fCrystalHalfLength + fCouplerHalfThickness * 2 + fPMTWindowHalfThickness);
    // const auto translation3 = G4ThreeVector(0., offsetlength, fCrystalHalfLength + fCouplerHalfThickness * 2 + fPMTWindowHalfThickness * 2 + fSiPMHalfThickness);

    const auto Transform =
        [&fCrystalLength, crystalTail = G4ThreeVector(0, 0, fCrystalLength / 2)](double transformDistance) {
            return G4Translate3D{crystalTail + G4ThreeVector(0, 0, transformDistance)};
        };

    // const G4int nsect = 5;
    // std::vector<G4TwoVector> polygon(nsect);
    // G4double ang = twopi / nsect;
    // for (int i = 0; i < nsect; ++i) {
    //     G4double phi = i * ang;
    //     G4double cosphi = std::cos(phi);
    //     G4double sinphi = std::sin(phi);
    //     polygon[i].set(fCrystalHalfWidth * cosphi, fCrystalHalfWidth * sinphi);
    // }

    // G4TwoVector offsetA(0, 0), offsetB(0, offsetlength);
    // G4double scaleA = 0.6, scaleB = 1;

    // const auto crystalSV = new G4ExtrudedSolid("Extruded", polygon, fCrystalHalfLength, offsetA, scaleA, offsetB, scaleB);

    const auto crystalSV = new G4Box("crystal", fCrystalWidth / 2, fCrystalWidth / 2, fCrystalLength / 2);

    //========================================== CsI(Tl) ============================================
    const auto crystalLV = new G4LogicalVolume{crystalSV, csI, "crystal"};
    //========================================== LaBr3(Ce) ==========================================
    // const auto crystalLV = new G4LogicalVolume{crystalSV, labr, "crystal"};
    //========================================== LYSO(Ce) ===========================================
    // const auto crystalLV = new G4LogicalVolume{crystalSV, lyso, "crystal"};
    //===============================================================================================
    const auto crystalPV = new G4PVPlacement(G4Transform3D{}, crystalLV, "crystal", worldLV, false, 0, true);

    // const auto couplerSV = new G4Tubs("coupler", 0, fCouplerWidth / 2, fCouplerThickness / 2, 0, 2 * pi);
    const auto headCouplerSV = new G4Box("coupler", fHeadCouplerWidth / 2, fHeadCouplerWidth / 2, fCouplerThickness / 2);
    const auto headCouplerLV = new G4LogicalVolume(headCouplerSV, siliconeOil, "CrystalCoupler");
    const auto headCouplerPV = new G4PVPlacement(Transform(fCouplerThickness / 2),
                                                 headCouplerLV, "CrystalCoupler", worldLV, false, 0, true);

    const auto guideSV = new G4Trd("guide", fLightGuideHeadWidth / 2, fLightGuideTailWidth / 2, fLightGuideHeadWidth / 2, fLightGuideTailWidth / 2, fLightGuideThickness / 2);
    const auto guideLV = new G4LogicalVolume(guideSV, acrylic, "guide");
    const auto guidePV = new G4PVPlacement(Transform(fCouplerThickness + fLightGuideThickness / 2),
                                           guideLV, "guide", worldLV, false, 0, true);

    const auto tailCouplerSV = new G4Box("coupler", fTailCouplerWidth / 2, fTailCouplerWidth / 2, fCouplerThickness / 2);
    const auto tailCouplerLV = new G4LogicalVolume(tailCouplerSV, siliconeOil, "CrystalCoupler");
    const auto tailouplerPV = new G4PVPlacement(Transform(fCouplerThickness + fLightGuideThickness + fCouplerThickness / 2),
                                                tailCouplerLV, "CrystalCoupler", worldLV, false, 0, true);

    // const auto windowSV = new G4Tubs("window", 0, fPMTWindowDiameter, fPMTWindowHalfThickness, 0, 2 * pi);
    // const auto windowLV = new G4LogicalVolume(windowSV, glass, "window");
    // const auto windowPV = new G4PVPlacement{nullptr, translation2, windowLV, "window", worldLV, false, 0, true};

    const auto sipmSV = new G4Box("sipm", fSiPMWidth / 2, fSiPMWidth / 2, fSiPMThickness / 2);
    const auto sipmLV = new G4LogicalVolume(sipmSV, silicon, "sipm");
    new G4PVPlacement(Transform(fCouplerThickness + fLightGuideThickness + fCouplerThickness + fSiPMThickness / 2),
                      sipmLV, "sipm", worldLV, false, 0, true);

    // Define Surface

    const auto rfSurface = new G4OpticalSurface("reflector", unified, polished, dielectric_metal);
    rfSurface->SetMaterialPropertiesTable(rfSurfacePropertiesTable);
    new G4LogicalBorderSurface("rfSurface", crystalPV, worldPV, rfSurface);
    new G4LogicalBorderSurface("rfSurface", guidePV, worldPV, rfSurface);

    const auto couplerSurface = new G4OpticalSurface("coupler", unified, polished, dielectric_dielectric);
    couplerSurface->SetMaterialPropertiesTable(couplerSurfacePropertiesTable);
    new G4LogicalBorderSurface("couplerSurface", crystalPV, headCouplerPV, couplerSurface);

    const auto airPaintSurface = new G4OpticalSurface("Paint", unified, polished, dielectric_dielectric);
    airPaintSurface->SetMaterialPropertiesTable(airPaintSurfacePropertiesTable);
    new G4LogicalBorderSurface("AirPaintSurface", worldPV, crystalPV, airPaintSurface);
    new G4LogicalBorderSurface("AirPaintSurface", worldPV, guidePV, airPaintSurface);

    const auto cathodeSurface = new G4OpticalSurface("Cathode", unified, polished, dielectric_metal);
    cathodeSurface->SetMaterialPropertiesTable(cathodeSurfacePropertiesTable);
    new G4LogicalSkinSurface("cathodeSkinSurface", sipmLV, cathodeSurface);

    // auto cathodeSurface = new G4OpticalSurface("Cathode", glisur, polished, dielectric_metal);
    // new G4LogicalSkinSurface("sipmSkinSurface", sipmLV, cathodeSurface);
    // cathodeSurfacePropertiesTable->AddProperty(
    //     "REFLECTIVITY",
    //     cathodeSurfacePropertiesEnergy,
    //     cathodeSurfaceProperties["REFLECTIVITY"]);
    // cathodeSurfacePropertiesTable->AddProperty(
    //     "EFFICIENCY",
    //     cathodeSurfacePropertiesEnergy,
    //     cathodeSurfacePropertiesEfficiency);
    // cathodeSurface->SetMaterialPropertiesTable(cathodeSurfacePropertiesTable);
    // cathodeSurfacePropertiesTable->DumpTable();

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