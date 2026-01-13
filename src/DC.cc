//-------- GEANT4
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4PVReplica.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4SDParticleFilter.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VPVParameterisation.hh"
//#include "G4PVParameterised.hh"
//#include "G4OpBoundaryProcess.hh"
#include "G4SDManager.hh"
//#include "G4ThreeVector.hh"
//#include "G4RotationMatrix.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"

#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4UnitsTable.hh"

#include "DC.hh"

DC::DC(double density, double lunghezza_collimatore)
 :  experimentalHall_log(0),GAS_log(0),PPAC_log(0),SiPM_log(0),collf_log(0),coll_log(0),cathode_log(0),anode_log(0),cathodeAl_log(0),anodeAl_log(0),
    experimentalHall_phys(0),GAS_phys(0),PPAC_phys(0),SiPM_phys(0),collf_phys(0),coll_phys(0),cathode_phys(0),anode_phys(0),cathodeAl_phys(0),anodeAl_phys(0),reflectMPT(nullptr), absorbMPT(nullptr), cf4SiSurface(nullptr)
{
  dens = density;
  collimatore = lunghezza_collimatore;
}

DC::~DC()
{


   // Clean up optical surfaces
    if (reflectMPT) delete reflectMPT;
    if (absorbMPT) delete absorbMPT; 
    if (cf4SiSurface) delete cf4SiSurface;
}

G4VPhysicalVolume* DC::Construct()
{
  DefineMaterials();
  ConstructLaboratory();

  return experimentalHall_phys;
}

void
DC::DefineMaterials()
{

  G4NistManager* pNistManager = G4NistManager::Instance();
  G4double a, z,density;
  G4int nelements;
  G4int ncomponents, natoms;


  //================================== elements ===================================
  //G4Element *H  = new G4Element("Hydrogen","H",1.,1.0079*g/mole);
  G4Element *C  = new G4Element("Carbon","C",6.,12.011*g/mole);
  G4Element *N  = new G4Element("Nitrogen","N",7.,14.007*g/mole);
  G4Element *O  = new G4Element("Oxygen","O",8.,15.999*g/mole);
  G4Element *F  = new G4Element("Fluorine","F",9.,18.998*g/mole);
  //G4Element* Ar = new G4Element("Argon","Ar",18.,39.95*g/mole);
  G4Element* Au = new G4Element("Gold","Au",79.,196.9665*g/mole);

  // boron
  G4Isotope* B10 = new G4Isotope("B10", 5, 10);
  G4Isotope* B11 = new G4Isotope("B11", 5, 11);
  G4Element* B = new G4Element("Boron", "B", ncomponents=2);
  B->AddIsotope(B10, 19.9*perCent);
  B->AddIsotope(B11, 80.1*perCent);
  G4Material* boron = new G4Material("boron", 2.46*g/cm3, ncomponents=1, kStateSolid,293*kelvin, 1*atmosphere);
  boron->AddElement(B, natoms=1);


  //----------------------------------- vacuum ------------------------------------
  //G4Material* Vacuum = pNistManager->FindOrBuildMaterial("G4_Galactic");
  G4Material *Vacuum = new G4Material("Vacuum", 1.e-20*g/cm3, 2, kStateGas);
  Vacuum->AddElement(N, 0.755);
  Vacuum->AddElement(O, 0.245);
  //gMan->SetMaterialAppearance(vacuum,-1);
  //::::::::::::::::::::::::::::::: Define the lead material:::::::::::::
  leadMaterial = new G4Material("Lead", 82, 207.2 * g/mole, 11.35 * g/cm3);

  //----------------------------------- Aluminium ------------------------------------
 
   
   Aluminium = new G4Material("Al", z = 13., a = 26.98 * g / mole,
                        density = 2.7 * g / cm3);
 

  //....=========================Polyethylene=======================


   // polyethilene
  G4Element* Hpe = new G4Element("TS_H_of_Polyethylene", "H", 1, 1.0079*g/mole);
  G4Element* Cpe = new G4Element("Carbon", "C", 6, 12.01*g/mole);
  polyethylene = new G4Material("polyethylene", 0.93*g/cm3, ncomponents=2, kStateSolid, 293*kelvin, 1*atmosphere);
  polyethylene->AddElement(Hpe, natoms=4);
  polyethylene->AddElement(Cpe, natoms=2);

  //================== borated polyethilene====================================
  b_polyethylene = new G4Material("b_polyethylene",0.94*g/cm3,ncomponents=4,kStateSolid,293*kelvin,1*atmosphere);
  b_polyethylene->AddElement(Hpe, 11.6*perCent);
  b_polyethylene->AddElement(Cpe, 61.2*perCent);
  b_polyethylene->AddElement(B, 5*perCent);
  b_polyethylene->AddElement(O, 22.2*perCent);
  

  //=================silicon material=============================================
   silicon = pNistManager->FindOrBuildMaterial("G4_Si");

  // standard density = 3.72 mg/cm3
  //dens= 0.046*atmosphere; //35Torr
  //----------------------------------- CarbonTetrafluoride ------------------------
  G4Material* CF4 = new G4Material("CF4", dens*mg/cm3,2,kStateGas);
  CF4->AddElement(C,1);
  CF4->AddElement(F,4);
  //gMan->SetMaterialAppearance(CF4,kBlue-10,0.1);

  const G4int iNbEntries = 300;
  std::vector<G4double>CF4PhotonMomentum  = {6.2*eV,6.138613861*eV,6.078431373*eV,6.019417476*eV,5.961538462*eV,5.904761905*eV,5.849056604*eV,5.794392523*eV,5.740740741*eV,5.688073394*eV,5.636363636*eV,5.585585586*eV,
						5.535714286*eV,5.486725664*eV,5.438596491*eV,5.391304348*eV,5.344827586*eV,5.299145299*eV,5.254237288*eV,5.210084034*eV,5.166666667*eV,5.123966942*eV,5.081967213*eV,5.040650407*eV,5*eV,4.96*eV,4.920634921*eV,
						4.881889764*eV,4.84375*eV,4.80620155*eV,4.769230769*eV,4.732824427*eV,4.696969697*eV,4.661654135*eV,4.626865672*eV,4.592592593*eV,4.558823529*eV,4.525547445*eV,4.492753623*eV,4.460431655*eV,
						4.428571429*eV,4.397163121*eV,4.366197183*eV,4.335664336*eV,4.305555556*eV,4.275862069*eV,4.246575342*eV,4.217687075*eV,4.189189189*eV,4.161073826*eV,4.133333333*eV,4.105960265*eV,4.078947368*eV,4.052287582*eV,
						4.025974026*eV,4*eV,3.974358974*eV,3.949044586*eV,3.924050633*eV,3.899371069*eV,3.875*eV,3.850931677*eV,3.827160494*eV,3.803680982*eV,3.780487805*eV,3.757575758*eV,3.734939759*eV,
						3.71257485*eV,3.69047619*eV,3.668639053*eV,3.647058824*eV,3.625730994*eV,3.604651163*eV,3.583815029*eV,3.563218391*eV,3.542857143*eV,3.522727273*eV,3.502824859*eV,3.483146067*eV,
						3.463687151*eV,3.444444444*eV,3.425414365*eV,3.406593407*eV,3.387978142*eV,3.369565217*eV,3.351351351*eV,3.333333333*eV,3.315508021*eV,3.29787234*eV,3.28042328*eV,3.263157895*eV,
						3.246073298*eV,3.229166667*eV,3.212435233*eV,3.195876289*eV,3.179487179*eV,3.163265306*eV,3.147208122*eV,3.131313131*eV,3.115577889*eV,3.1*eV,3.084577114*eV,3.069306931*eV,3.054187192*eV,
						3.039215686*eV,3.024390244*eV,3.009708738*eV,2.995169082*eV,2.980769231*eV,2.966507177*eV,2.952380952*eV,2.938388626*eV,2.924528302*eV,2.910798122*eV,2.897196262*eV,2.88372093*eV,
						2.87037037*eV,2.857142857*eV,2.844036697*eV,2.831050228*eV,2.818181818*eV,2.805429864*eV,2.792792793*eV,2.780269058*eV,2.767857143*eV,2.755555556*eV,2.743362832*eV,
						2.731277533*eV,2.719298246*eV,2.707423581*eV,2.695652174*eV,2.683982684*eV,2.672413793*eV,2.660944206*eV,2.64957265*eV,2.638297872*eV,2.627118644*eV,2.616033755*eV,2.605042017*eV,
						2.594142259*eV,2.583333333*eV,2.572614108*eV,2.561983471*eV,2.551440329*eV,2.540983607*eV,2.530612245*eV,2.520325203*eV,2.510121457*eV,2.5*eV,2.489959839*eV,2.48*eV,2.470119522*eV,2.46031746*eV,
						2.450592885*eV,2.440944882*eV,2.431372549*eV,2.421875*eV,2.412451362*eV,2.403100775*eV,2.393822394*eV,2.384615385*eV,2.375478927*eV,2.366412214*eV,2.357414449*eV,
						2.348484848*eV,2.339622642*eV,2.330827068*eV,2.322097378*eV,2.313432836*eV,2.304832714*eV,2.296296296*eV,2.287822878*eV,2.279411765*eV,2.271062271*eV,2.262773723*eV,2.254545455*eV,2.246376812*eV,2.238267148*eV,
						2.230215827*eV,2.222222222*eV,2.214285714*eV,2.206405694*eV,2.19858156*eV,2.190812721*eV,2.183098592*eV,2.175438596*eV,2.167832168*eV,2.160278746*eV,2.152777778*eV,2.14532872*eV,2.137931034*eV,2.130584192*eV,
						2.123287671*eV,2.116040956*eV,2.108843537*eV,2.101694915*eV,2.094594595*eV,2.087542088*eV,2.080536913*eV,2.073578595*eV,2.066666667*eV,2.059800664*eV,2.052980132*eV,2.04620462*eV,2.039473684*eV,2.032786885*eV,
						2.026143791*eV,2.019543974*eV,2.012987013*eV,2.006472492*eV,2*eV,1.993569132*eV,1.987179487*eV,1.980830671*eV,1.974522293*eV,1.968253968*eV,1.962025316*eV,1.955835962*eV,
						1.949685535*eV,1.943573668*eV,1.9375*eV,1.931464174*eV,1.925465839*eV,1.919504644*eV,1.913580247*eV,1.907692308*eV,1.901840491*eV,1.896024465*eV,1.890243902*eV,1.88449848*eV,1.878787879*eV,1.873111782*eV,
						1.86746988*eV,1.861861862*eV,1.856287425*eV,1.850746269*eV,1.845238095*eV,1.839762611*eV,1.834319527*eV,1.828908555*eV,1.823529412*eV,1.818181818*eV,
						1.812865497*eV,1.807580175*eV,1.802325581*eV,1.797101449*eV,1.791907514*eV,1.786743516*eV,1.781609195*eV,1.776504298*eV,1.771428571*eV,1.766381766*eV,1.761363636*eV,1.756373938*eV,1.751412429*eV,
						1.746478873*eV,1.741573034*eV,1.736694678*eV,1.731843575*eV,1.727019499*eV,1.722222222*eV,1.717451524*eV,1.712707182*eV,1.707988981*eV,1.703296703*eV,1.698630137*eV,1.693989071*eV,1.689373297*eV,1.684782609*eV,
						1.680216802*eV,1.675675676*eV,1.67115903*eV,1.666666667*eV,1.662198391*eV,1.657754011*eV,1.653333333*eV,1.64893617*eV,1.644562334*eV,1.64021164*eV,1.635883905*eV,1.631578947*eV,1.627296588*eV,1.623036649*eV,
						1.618798956*eV,1.614583333*eV,1.61038961*eV,1.606217617*eV,1.602067183*eV,1.597938144*eV,1.593830334*eV,1.58974359*eV,1.585677749*eV,1.581632653*eV,1.577608142*eV,1.573604061*eV,1.569620253*eV,
						1.565656566*eV,1.561712846*eV,1.557788945*eV,1.553884712*eV};
    // Sort the array in ascending order
  std::sort(CF4PhotonMomentum.begin(), CF4PhotonMomentum.end());

  std::vector<G4double>CF4Scintillation_Fast     = {0.0029,0.0029,0.0017,0.0024,0.0018,0.0011,0.0027,0.0009,0.0003,0.0019,0.0030,0.0024,0.0023,0.0036,0.0039,0.0056,
						0.0049,0.0061,0.0053,0.0052,0.0056,0.0064,0.0072,0.0064,0.0080,0.0071,0.0056,0.0069,0.0053,0.0070,0.0060,0.0057,0.0071,0.0066,0.0066,
						0.0055,0.0082,0.0076,0.0093,0.0089,0.0106,0.0109,0.0105,0.0102,0.0120,0.0121,0.0102,0.0097,0.0120,0.0126,0.0097,0.0103,0.0097,0.0084,
						0.0119,0.0112,0.0096,0.0171,0.0235,0.0078,0.0089,0.0071,0.0065,0.0074,0.0073,0.0074,0.0074,0.0080,0.0143,0.0522,0.0069,0.0076,0.0042,
						0.0059,0.0039,0.0053,0.0054,0.0185,0.0077,0.0599,0.0048,0.0034,0.0041,0.0041,0.0047,0.0059,0.0046,0.0065,0.0128,0.0037,0.0167,0.0053,
						0.0038,0.0042,0.0046,0.0032,0.0037,0.0073,0.0049,0.0067,0.0116,0.0054,0.0077,0.0111,0.0042,0.0043,0.0037,0.0046,0.0041,0.0028,0.0055,
						0.0031,0.0048,0.0057,0.0056,0.0035,0.0039,0.0068,0.0051,0.0037,0.0054,0.0048,0.0061,0.0033,0.0050,0.0052,0.0047,0.0014,0.0043,0.0041,
						0.0023,0.0062,0.0036,0.0038,0.0039,0.0043,0.0049,0.0049,0.0036,0.0048,0.0039,0.0023,0.0035,0.0025,0.0036,0.0010,0.0044,0.0013,0.0041,
						0.0021,0.0016,0.0046,0.0040,0.0034,0.0027,0.0026,0.0034,0.0004,0.0037,0.0004,0.0036,0.0029,0.0029,0.0036,0.0055,0.0034,0.0034,0.0025,
						0.0028,0.0055,0.0064,0.0037,0.0029,0.0047,0.0058,0.0040,0.0062,0.0055,0.0029,0.0067,0.0070,0.0080,0.0060,0.0094,0.0082,0.0072,0.0089,
						0.0117,0.0102,0.0134,0.0131,0.0131,0.0120,0.0135,0.0096,0.0107,0.0179,0.0210,0.0172,0.0165,0.0167,0.0176,0.0137,0.0196,0.0217,0.0175,
						0.0223,0.0192,0.0222,0.0188,0.0184,0.0183,0.0156,0.0098,0.0198,0.0268,0.0188,0.0236,0.0208,0.0171,0.0229,0.0228,0.0227,0.0204,0.0184,
						0.0190,0.0185,0.0145,0.0138,0.0122,0.0180,0.0132,0.0146,0.0087,0.0039,0.0147,0.0000,0.0000,0.0137,0.0084,0.0094,0.0114,0.0078,0.0100,
						0.0069,0.0055,0.0164,0.0113,0.0148,0.0053,0.0054,0.0065,0.0092,0.0000,0.0047,0.0000,0.0071,0.0000,0.0057,0.0063,0.0064,0.0050,0.0077,
						0.0034,0.0025,0.0000,0.0041,0.0025,0.0019,0.0042,0.0030,0.0000,0.0030,0.0000,0.0000,0.0000,0.0027,0.0000,0.0000,0.0000,0.0000,0.0006,
						0.0051,0.0083,0.0000,0.0000,0.0064,0.0003,0.0002,0.0074,0.0038,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000};
   std::vector<G4double>CF4Scintillation_Slow   = {0.0029,0.0029,0.0017,0.0024,0.0018,0.0011,0.0027,0.0009,0.0003,0.0019,0.0030,0.0024,0.0023,0.0036,0.0039,0.0056,
						0.0049,0.0061,0.0053,0.0052,0.0056,0.0064,0.0072,0.0064,0.0080,0.0071,0.0056,0.0069,0.0053,0.0070,0.0060,0.0057,0.0071,0.0066,0.0066,
						0.0055,0.0082,0.0076,0.0093,0.0089,0.0106,0.0109,0.0105,0.0102,0.0120,0.0121,0.0102,0.0097,0.0120,0.0126,0.0097,0.0103,0.0097,0.0084,
						0.0119,0.0112,0.0096,0.0171,0.0235,0.0078,0.0089,0.0071,0.0065,0.0074,0.0073,0.0074,0.0074,0.0080,0.0143,0.0522,0.0069,0.0076,0.0042,
						0.0059,0.0039,0.0053,0.0054,0.0185,0.0077,0.0599,0.0048,0.0034,0.0041,0.0041,0.0047,0.0059,0.0046,0.0065,0.0128,0.0037,0.0167,0.0053,
						0.0038,0.0042,0.0046,0.0032,0.0037,0.0073,0.0049,0.0067,0.0116,0.0054,0.0077,0.0111,0.0042,0.0043,0.0037,0.0046,0.0041,0.0028,0.0055,
						0.0031,0.0048,0.0057,0.0056,0.0035,0.0039,0.0068,0.0051,0.0037,0.0054,0.0048,0.0061,0.0033,0.0050,0.0052,0.0047,0.0014,0.0043,0.0041,
						0.0023,0.0062,0.0036,0.0038,0.0039,0.0043,0.0049,0.0049,0.0036,0.0048,0.0039,0.0023,0.0035,0.0025,0.0036,0.0010,0.0044,0.0013,0.0041,
						0.0021,0.0016,0.0046,0.0040,0.0034,0.0027,0.0026,0.0034,0.0004,0.0037,0.0004,0.0036,0.0029,0.0029,0.0036,0.0055,0.0034,0.0034,0.0025,
						0.0028,0.0055,0.0064,0.0037,0.0029,0.0047,0.0058,0.0040,0.0062,0.0055,0.0029,0.0067,0.0070,0.0080,0.0060,0.0094,0.0082,0.0072,0.0089,
						0.0117,0.0102,0.0134,0.0131,0.0131,0.0120,0.0135,0.0096,0.0107,0.0179,0.0210,0.0172,0.0165,0.0167,0.0176,0.0137,0.0196,0.0217,0.0175,
						0.0223,0.0192,0.0222,0.0188,0.0184,0.0183,0.0156,0.0098,0.0198,0.0268,0.0188,0.0236,0.0208,0.0171,0.0229,0.0228,0.0227,0.0204,0.0184,
						0.0190,0.0185,0.0145,0.0138,0.0122,0.0180,0.0132,0.0146,0.0087,0.0039,0.0147,0.0000,0.0000,0.0137,0.0084,0.0094,0.0114,0.0078,0.0100,
						0.0069,0.0055,0.0164,0.0113,0.0148,0.0053,0.0054,0.0065,0.0092,0.0000,0.0047,0.0000,0.0071,0.0000,0.0057,0.0063,0.0064,0.0050,0.0077,
						0.0034,0.0025,0.0000,0.0041,0.0025,0.0019,0.0042,0.0030,0.0000,0.0030,0.0000,0.0000,0.0000,0.0027,0.0000,0.0000,0.0000,0.0000,0.0006,
						0.0051,0.0083,0.0000,0.0000,0.0064,0.0003,0.0002,0.0074,0.0038,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000};


  const G4int iNbEntries_1 = 3;
  G4double CF4PhotonMomentum_1[iNbEntries_1] = {200*eV,500*eV,798*eV};
  G4double CF4RefractiveIndex[iNbEntries_1]  = {1.004,1.004,1.004};
  G4double CF4AbsorbtionLength[iNbEntries_1] = {100.*cm, 100.*cm, 100.*cm};
  G4double CF4ScatteringLength[iNbEntries_1] = {30.*cm,  30.*cm,  30.*cm};
  G4MaterialPropertiesTable *CF4PropertiesTable = new G4MaterialPropertiesTable();
  CF4PropertiesTable->AddProperty("SCINTILLATIONCOMPONENT1", CF4PhotonMomentum, CF4Scintillation_Fast, iNbEntries);
  CF4PropertiesTable->AddProperty("SCINTILLATIONCOMPONENT2", CF4PhotonMomentum, CF4Scintillation_Slow, iNbEntries);
  CF4PropertiesTable->AddProperty("RINDEX", CF4PhotonMomentum_1, CF4RefractiveIndex, iNbEntries_1);
  CF4PropertiesTable->AddProperty("ABSLENGTH", CF4PhotonMomentum_1, CF4AbsorbtionLength, iNbEntries_1);
  CF4PropertiesTable->AddProperty("RAYLEIGH", CF4PhotonMomentum_1, CF4ScatteringLength, iNbEntries_1);
  CF4PropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 2500./keV,true);  // for electron recoil
  CF4PropertiesTable->AddConstProperty("RESOLUTIONSCALE", 1.0);
  CF4PropertiesTable->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 3.*ns,true);
  CF4PropertiesTable->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 10.*ns,true);
  CF4PropertiesTable->AddConstProperty("YIELDRATIO", 1.0,true);
  CF4->SetMaterialPropertiesTable(CF4PropertiesTable);

  //------------------------------------ teflon -----------------------------------

  const G4int TNbEntries = 4;
  G4Material* Teflon = new G4Material("Teflon", 2.2*g/cm3, 2, kStateSolid);
  Teflon->AddElement(C, 0.240183);
  Teflon->AddElement(F, 0.759817);

  G4double pdTeflonPhotonMomentum[TNbEntries]  = {1.0*eV,6.91*eV, 6.98*eV, 7.05*eV};
  G4double pdTeflonRefractiveIndex[TNbEntries] = {1.34,1.34,1.34,1.34};
  G4double pdTeflonReflectivity[TNbEntries]    = {0.9,0.9,0.9,0.9};
  G4double pdTeflonAbsLength[TNbEntries]       = {0.2*mm,0.2*mm,0.2*mm,0.2*mm};
  G4double pdTeflonScatteringLength[TNbEntries] = {30*cm,30*cm,30*cm,30*cm};

  G4MaterialPropertiesTable *pTeflonPropertiesTable = new G4MaterialPropertiesTable();

  pTeflonPropertiesTable->AddProperty("RINDEX", pdTeflonPhotonMomentum, pdTeflonRefractiveIndex, TNbEntries);

  pTeflonPropertiesTable->AddProperty("ABSLENGTH", pdTeflonPhotonMomentum, pdTeflonAbsLength, TNbEntries);  
  pTeflonPropertiesTable->AddProperty("RAYLEIGH", pdTeflonPhotonMomentum, pdTeflonScatteringLength, TNbEntries);
  pTeflonPropertiesTable->AddProperty("REFLECTIVITY", pdTeflonPhotonMomentum, pdTeflonReflectivity, TNbEntries);
  Teflon->SetMaterialPropertiesTable(pTeflonPropertiesTable);

// 
// ------------------------------------ Polypropilene ------------------------------------

pNistManager->FindOrBuildMaterial("G4_POLYPROPYLENE");

// ------------------------------------ Gold - Au ------------------------------------
//
const G4int iNbEntries_Au = 4;
G4Material* Gold = new G4Material("Gold",19.30*g/cm3,1,kStateSolid);
Gold->AddElement(Au, 1.0);

G4double AuPM[iNbEntries_Au]  = {1*eV,  2.5001*eV, 2.5001*eV, 7*eV};
G4double AuReflectivity[iNbEntries_Au] = {0.3, 0.3, 1, 1};
G4double AuRefractiveIndex[iNbEntries_Au] = {1.25,0.86,0.75,3.56};
G4double AuAbsLength[iNbEntries_Au] = {1e-5*mm,1e-5*mm,1e-5*mm,1e-5*mm};
G4double AuScatteringLength[iNbEntries_Au] = {30*cm,30*cm,30*cm,30*cm};

G4MaterialPropertiesTable *AuPropertiesTable = new G4MaterialPropertiesTable();

AuPropertiesTable->AddProperty("RINDEX", AuPM, AuRefractiveIndex, iNbEntries_Au);
AuPropertiesTable->AddProperty("ABSLENGTH", AuPM, AuAbsLength, iNbEntries_Au);  
AuPropertiesTable->AddProperty("RAYLEIGH", AuPM, AuScatteringLength, iNbEntries_Au);
AuPropertiesTable->AddProperty("REFLECTIVITY", AuPM, AuReflectivity, iNbEntries_Au);
////::::::::::::::::::: SURFACE PROPERTIES:::::::::::::::::::::
//REFLECTING SURFACE (Aluminium cathode/anode)
reflectMPT = new G4MaterialPropertiesTable();

// Aluminium reflectivity - high in visible/UV range
const G4int numEntriesAl = 11;
G4double alPhotonEnergy[numEntriesAl] = {1.13*eV, 1.24*eV, 1.38*eV, 1.55*eV, 1.77*eV, 2.07*eV, 
                                       2.48*eV, 3.10*eV, 4.13*eV, 6.20*eV, 6.53*eV};

// Aluminium reflectivity (typical values - high reflectivity)
G4double alReflectivity[numEntriesAl] = {
    0.92, 0.92, 0.91, 0.90, 0.89, 0.87,
    0.85, 0.82, 0.78, 0.70, 0.65
};

// Zero efficiency for metals (no transmission)
G4double alEfficiency[numEntriesAl] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

reflectMPT->AddProperty("REFLECTIVITY", alPhotonEnergy, alReflectivity, numEntriesAl);
reflectMPT->AddProperty("EFFICIENCY", alPhotonEnergy, alEfficiency, numEntriesAl);


//ABSORBING SURFACE (Teflon collimators)
absorbMPT = new G4MaterialPropertiesTable();

// Teflon properties - low reflectivity, acts as absorber
G4double teflonReflectivity[numEntriesAl] = {
    0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
    0.05, 0.05, 0.05, 0.05, 0.05  // 5% reflectivity
};

// Zero efficiency - absorbs all transmitted light
G4double teflonEfficiency[numEntriesAl] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

absorbMPT->AddProperty("REFLECTIVITY", alPhotonEnergy, teflonReflectivity, numEntriesAl);
absorbMPT->AddProperty("EFFICIENCY", alPhotonEnergy, teflonEfficiency, numEntriesAl);

//:::::::::::::::::::::::::::optical properties for silicon::::::::::::::::::::::::::::::::::::::::::::::
 // =========================================================================
    // BULK MATERIAL PROPERTIES TABLE
    // =========================================================================
    
    G4MaterialPropertiesTable* siMPT = new G4MaterialPropertiesTable();
    const G4int iNbEntries_Si = 11;
    
    G4double SiPM[iNbEntries_Si]  = {1.13*eV, 1.24*eV, 1.38*eV, 1.55*eV, 1.77*eV, 2.07*eV, 2.48*eV, 3.10*eV, 4.13*eV, 6.20*eV, 6.53*eV};// Define all photon energies (in eV)
    G4double SiRefractiveIndex[iNbEntries_Si] = {3.60, 3.64, 3.68, 3.73, 3.81, 3.94, 4.30, 5.57, 4.63, 1.78, 1.65};     // Refractive index (n)
    //G4double SiRefractiveIndex[iNbEntries_Si] = {1.25,0.86,0.75,3.56};
    G4double SiAbsLength[iNbEntries_Si] = {10.0000*mm, 5.00000*mm, 2.00000*mm, 0.50000*mm, 0.10000*mm, 0.03000*mm, 0.01000*mm, 0.00100*mm, 0.00010*mm, 0.00001*mm, 0.00001*mm};       // Absorption length
    G4double SiScatteringLength[iNbEntries_Si] = {10000.0*m, 10000.0*m, 10000.0*m, 10000.0*m, 10000.0*m, 1000.0*m, 100.0*m, 10.0*m, 1.0*m, 0.1*m, 0.1*m};  // Rayleigh scattering length (typically very long for pure silicon)


    siMPT->AddProperty("RINDEX", SiPM, SiRefractiveIndex, iNbEntries_Si);
    siMPT->AddProperty("ABSLENGTH", SiPM, SiAbsLength, iNbEntries_Si);   
    siMPT->AddProperty("RAYLEIGH", SiPM, SiScatteringLength, iNbEntries_Si);
    //SiMPT->AddProperty("REFLECTIVITY", SiPM, SiReflectivity, iNbEntries_Si);
    siMPT->AddConstProperty("RESOLUTIONSCALE", 1.0);

    silicon->SetMaterialPropertiesTable(siMPT);

    // =========================================================================
        // CF4 - SILICON OPTICAL BOUNDARY SURFACE
    // =========================================================================

       // Create optical surface for CF4-Si interface
    cf4SiSurface = new G4OpticalSurface("CF4_Si_Surface");
    cf4SiSurface->SetType(dielectric_dielectric);
    cf4SiSurface->SetFinish(polished);
    cf4SiSurface->SetModel(glisur);

    // CF4 refractive index (approximately 1.0005 - close to vacuum)
    G4double cf4Rindex[iNbEntries_Si] = {1.0005, 1.0005, 1.0005, 1.0005, 1.0005, 1.0005,1.0005, 1.0005, 1.0005, 1.0005, 1.0005};
    // Calculate reflectivity using Fresnel equations for CF4-Si interface
    G4double reflectivity[iNbEntries_Si];
    for (int i=0; i<iNbEntries_Si; i++) {
        G4double n1 = cf4Rindex[i]; // CF4 refractive index
        G4double n2 = SiRefractiveIndex[i];  // Silicon refractive index
        // Fresnel equation for normal incidence: R = [(n2-n1)/(n2+n1)]^2
        reflectivity[i] = ((n2 - n1) / (n2 + n1)) * ((n2 - n1) / (n2 + n1));
    }


    // Typical silicon photodetector efficiency 
    G4double efficiency[iNbEntries_Si] = {
     0.10,  // 1100 nm - lower efficiency in IR
     0.15,  // 1000 nm
     0.20,  // 900 nm
     0.25,  // 800 nm  
     0.30,  // 700 nm
     0.40,  // 600 nm
     0.60,  // 500 nm - peak efficiency in green
     0.80,  // 400 nm - high efficiency in blue/UV
     0.70,  // 300 nm - decreasing in deep UV
     0.50,  // 200 nm
     0.30   // 190 nm
    };

    G4MaterialPropertiesTable* cf4SiSurfMPT = new G4MaterialPropertiesTable();
    cf4SiSurfMPT->AddProperty("REFLECTIVITY", SiPM, reflectivity, iNbEntries_Si);
    cf4SiSurfMPT->AddProperty("EFFICIENCY", SiPM, efficiency, iNbEntries_Si);

    cf4SiSurface->SetMaterialPropertiesTable(cf4SiSurfMPT);

}


void DC::Sphereball( G4double position) {
     G4double minSphereradius = 0*cm;
     G4double maxSphereradius = 20.0*cm;
     G4Sphere* sphereball = new G4Sphere("sphereball", minSphereradius/2, maxSphereradius/2 , 0*deg,360*deg,0*deg,180*deg);
     G4LogicalVolume* sphereVolume = new G4LogicalVolume(sphereball, polyethylene, "Sphere");
     G4PVPlacement* spherePlacement  = new G4PVPlacement(0,
                                               G4ThreeVector(0.*mm, 0.*mm, position),
                                               sphereVolume,
                                               "Sphere",
                                               experimentalHall_log,
                                               false,
                                               0,true);
     G4VisAttributes* blue = new G4VisAttributes(G4Colour::Blue());

     blue->SetVisibility(true);
     blue->SetForceAuxEdgeVisible(true);


     sphereVolume->SetVisAttributes(blue);


}

void DC::ConstructValve(G4double position, G4double ballOffset) {
    G4NistManager* nist = G4NistManager::Instance();

    // Define dimensions
    G4double valveLength = 20.0 * cm;
    G4double valveRadius = 5.0 * cm;
    G4double passageRadius = 3.0 * cm;
    G4double ballRadius = 4.0 * cm;

    // Materials
    G4Material* steel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    G4Material* water = nist->FindOrBuildMaterial("G4_WATER");

    // **Step 1: Main Valve Body (Cylinder)**
    G4Tubs* valveBody = new G4Tubs("ValveBody", 0, valveRadius, valveLength/2, 0, 360*deg);

    // **Step 2: Hollow Passage (Subtracted Cylinder)**
    G4Tubs* flowPassage = new G4Tubs("FlowPassage", 0, passageRadius, valveLength, 0, 360*deg);
    G4SubtractionSolid* hollowValve = new G4SubtractionSolid("HollowValve", valveBody, flowPassage);

    // **Step 3: Ball Inside the Valve**
    G4Sphere* ball = new G4Sphere("ValveBall", 0, ballRadius, 0, 360*deg, 0, 180*deg);

    // Position the ball inside the valve (use ballOffset to "open/close")
    G4ThreeVector ballPosition(0, 0, ballOffset); // Adjust position dynamically
    G4SubtractionSolid* finalValve = new G4SubtractionSolid("FinalValve", hollowValve, ball, 0, ballPosition);

    // **Step 4: Create Logical Volume**
    G4LogicalVolume* valveLogic = new G4LogicalVolume(finalValve, steel, "ValveLogic");

    // **Step 5: Place in World**
    new G4PVPlacement(0, G4ThreeVector(0, 0, position), valveLogic, "ValvePhys", experimentalHall_log, false, 0, true);

/*    // **Step 6: Define the Fluid (Water) inside the passage**
    G4LogicalVolume* waterLogic = new G4LogicalVolume(flowPassage, water, "WaterLogic");
    new G4PVPlacement(0, G4ThreeVector(0, 0, position), waterLogic, "WaterPhys", experimentalHall_log, false, 0, true);
*/
     G4VisAttributes* blue = new G4VisAttributes(G4Colour::Blue());

     blue->SetVisibility(true);
     blue->SetForceAuxEdgeVisible(true);


     valveLogic->SetVisAttributes(blue);

     G4VisAttributes* red = new G4VisAttributes(G4Colour::Red());

     red->SetVisibility(true);
     red->SetForceAuxEdgeVisible(true);

/*
     waterLogic->SetVisAttributes(red);

*/
}



void DC::ConstructLaboratory()
{
  G4Material *Vacuum = G4Material::GetMaterial("Vacuum");
  G4Material *CF4 = G4Material::GetMaterial("CF4");
  G4Material *Teflon = G4Material::GetMaterial("Teflon");
  G4Material *PP = G4Material::GetMaterial("G4_POLYPROPYLENE");
  //G4Material *Silicon = G4Material::GetMaterial("Aluminium");

  //	========================================================================================================================
  //======================================================================================================================
  // Experimental hall (world volume) ==================================================================================

  G4double expHall_x = 5.0*m;
  G4double expHall_y = 5.0*m;
  G4double expHall_z = 5.0*m;
  G4Box* experimentalHall_box = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);
  experimentalHall_log = new G4LogicalVolume(experimentalHall_box,Vacuum,"expHall_log");
  experimentalHall_phys = new G4PVPlacement(0,G4ThreeVector(),experimentalHall_log,"expHall",0,false,0);
  experimentalHall_log->SetVisAttributes(G4VisAttributes::GetInvisible());

  //==========================MODERATOR AND COLLIMATOR===================================

  //The HDPE_block1
/*
  fblockSize = 10*cm;


  HDPE_Box1 = new G4Box("HDPE1",                             //its name
                   10*cm/2,10*cm/2,1*cm/2);   //its dimensions

  HDPE_LV1 = new G4LogicalVolume(HDPE_Box1,                     //its shape
                              polyethylene,                      //its material
                             "HDPE1");                  //its name

  HDPE_PV1 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,0,0.5*cm),            //at (0,0,0)
                             HDPE_LV1,                      //its logical volume
                            "HDPE1",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  //The HDPE_block2


  HDPE_Box2 = new G4Box("HDPE2",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV2 = new G4LogicalVolume(HDPE_Box2,                     //its shape
                              polyethylene,                      //its material
                             "HDPE2");                  //its name

  HDPE_PV2 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,10*cm,5*cm),            //at (0,0,0)
                             HDPE_LV2,                      //its logical volume
                            "HDPE2",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  //The HDPE_block3


  HDPE_Box3 = new G4Box("HDPE3",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV3 = new G4LogicalVolume(HDPE_Box3,                     //its shape
                              polyethylene,                      //its material
                             "HDPE3");                  //its name

  HDPE_PV3 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,-10*cm,5*cm),            //at (0,0,0)
                             HDPE_LV3,                      //its logical volume
                            "HDPE3",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  //The HDPE_block4


  HDPE_Box4 = new G4Box("HDPE4",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV4 = new G4LogicalVolume(HDPE_Box4,                     //its shape
                              polyethylene,                      //its material
                             "HDPE4");                  //its name

  HDPE_PV4 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,10*cm,15*cm),            //at (0,0,0)
                             HDPE_LV4,                      //its logical volume
                            "HDPE4",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  //The HDPE_block5


  HDPE_Box5 = new G4Box("HDPE5",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV5 = new G4LogicalVolume(HDPE_Box5,                     //its shape
                              polyethylene,                      //its material
                             "HDPE5");                  //its name

  HDPE_PV5 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,-10*cm,15*cm),            //at (0,0,0)
                             HDPE_LV5,                      //its logical volume
                            "HDPE5",                    //its name
                             experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number
 //The HDPE_block6


  HDPE_Box6 = new G4Box("HDPE6",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV6 = new G4LogicalVolume(HDPE_Box6,                     //its shape
                              polyethylene,                      //its material
                             "HDPE6");                  //its name

  HDPE_PV6 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,0,5*cm),            //at (0,0,0)
                             HDPE_LV6,                      //its logical volume
                            "HDPE6",                    //its name
                             experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number
  
   //The HDPE_block7


  HDPE_Box7 = new G4Box("HDPE7",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV7 = new G4LogicalVolume(HDPE_Box7,                     //its shape
                              polyethylene,                      //its material
                             "HDPE7");                  //its name


  HDPE_PV7 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,0,5*cm),            //at (0,0,0)
                             HDPE_LV7,                      //its logical volume
                            "HDPE7",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


//The HDPE_block8


  HDPE_Box8 = new G4Box("HDPE8",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV8 = new G4LogicalVolume(HDPE_Box8,                     //its shape
                              polyethylene,                      //its material
                             "HDPE8");                  //its name

  HDPE_PV8 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,10*cm,5*cm),            //at (0,0,0)
                             HDPE_LV8,                      //its logical volume
                            "HDPE8",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  //The HDPE_block9


  HDPE_Box9 = new G4Box("HDPE9",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV9 = new G4LogicalVolume(HDPE_Box9,                     //its shape
                              polyethylene,                      //its material
                             "HDPE9");                  //its name

  HDPE_PV9 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,-10*cm,5*cm),            //at (0,0,0)
                             HDPE_LV9,                      //its logical volume
                            "HDPE9",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

//The HDPE_block10

  HDPE_Box10 = new G4Box("HDPE10",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV10 = new G4LogicalVolume(HDPE_Box10,                     //its shape
                              polyethylene,                      //its material
                             "HDPE10");                  //its name

  HDPE_PV10 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,10*cm,5*cm),            //at (0,0,0)
                             HDPE_LV10,                      //its logical volume
                            "HDPE10",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

   //The HDPE_block11

  HDPE_Box11 = new G4Box("HDPE11",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV11 = new G4LogicalVolume(HDPE_Box11,                     //its shape
                              polyethylene,                      //its material
                             "HDPE11");                  //its name

  HDPE_PV11 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,-10*cm,5*cm),            //at (0,0,0)
                             HDPE_LV11,                      //its logical volume
                            "HDPE11",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


//The HDPE_block12


  HDPE_Box12 = new G4Box("HDPE12",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV12 = new G4LogicalVolume(HDPE_Box12,                     //its shape
                              polyethylene,                      //its material
                             "HDPE12");                  //its name

  HDPE_PV12 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,-10*cm,15*cm),            //at (0,0,0)
                             HDPE_LV12,                      //its logical volume
                            "HDPE12",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  //The HDPE_block13


  HDPE_Box13 = new G4Box("HDPE13",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV13 = new G4LogicalVolume(HDPE_Box13,                     //its shape
                              polyethylene,                      //its material
                             "HDPE13");                  //its name

  HDPE_PV13 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,10*cm,15*cm),            //at (0,0,0)
                             HDPE_LV13,                      //its logical volume
                            "HDPE13",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

//The HDPE_block14


  HDPE_Box14 = new G4Box("HDPE14",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV14 = new G4LogicalVolume(HDPE_Box14,                     //its shape
                              polyethylene,                      //its material
                             "HDPE14");                  //its name

  HDPE_PV14 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,-10*cm,15*cm),            //at (0,0,0)
                             HDPE_LV14,                      //its logical volume
                            "HDPE14",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  
  //The HDPE_block15

  HDPE_Box15 = new G4Box("HDPE15",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV15 = new G4LogicalVolume(HDPE_Box15,                     //its shape
                              polyethylene,                      //its material
                             "HDPE15");                  //its name

  HDPE_PV15 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,10*cm,15*cm),            //at (0,0,0)
                             HDPE_LV15,                      //its logical volume
                            "HDPE15",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

//The HDPE_block16


  HDPE_Box16 = new G4Box("HDPE16",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV16 = new G4LogicalVolume(HDPE_Box16,                     //its shape
                              polyethylene,                      //its material
                             "HDPE16");                  //its name

  HDPE_PV16 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,0*cm,15*cm),            //at (0,0,0)
                             HDPE_LV16,                      //its logical volume
                            "HDPE16",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  //The HDPE_block17

  HDPE_Box17 = new G4Box("HDPE17",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV17 = new G4LogicalVolume(HDPE_Box17,                     //its shape
                              polyethylene,                      //its material
                             "HDPE17");                  //its name

  HDPE_PV17 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,0*cm,15*cm),            //at (0,0,0)
                             HDPE_LV17,                      //its logical volume
                            "HDPE17",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

    //The HDPE_block18


  HDPE_Box18 = new G4Box("HDPE18",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  HDPE_LV18 = new G4LogicalVolume(HDPE_Box18,                     //its shape
                              polyethylene,                      //its material
                             "HDPE18");                  //its name

  HDPE_PV18 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,-10*cm,25*cm),            //at (0,0,0)
                             HDPE_LV18,                      //its logical volume
                            "HDPE18",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number
   //The HDPE_block19


  HDPE_Box19 = new G4Box("HDPE19",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  HDPE_LV19 = new G4LogicalVolume(HDPE_Box19,                     //its shape
                              polyethylene,                      //its material
                             "HDPE19");                  //its name

  HDPE_PV19 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,-10*cm,25*cm),            //at (0,0,0)
                             HDPE_LV19,                      //its logical volume
                            "HDPE19",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);

     //The HDPE_block20


  HDPE_Box20 = new G4Box("HDPE20",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  HDPE_LV20 = new G4LogicalVolume(HDPE_Box20,                     //its shape
                              polyethylene,                      //its material
                             "HDPE20");                  //its name

  HDPE_PV20 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,10*cm,25*cm),            //at (0,0,0)
                             HDPE_LV20,                      //its logical volume
                            "HDPE20",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number
   //The HDPE_block21

  HDPE_Box21 = new G4Box("HDPE21",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  HDPE_LV21 = new G4LogicalVolume(HDPE_Box21,                     //its shape
                              polyethylene,                      //its material
                             "HDPE21");                  //its name

  HDPE_PV21 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,10*cm,25*cm),            //at (0,0,0)
                             HDPE_LV21,                      //its logical volume
                            "HDPE21",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

//The HDPE_block22


  HDPE_Box22 = new G4Box("HDPE22",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  HDPE_LV22 = new G4LogicalVolume(HDPE_Box22,                     //its shape
                              polyethylene,                      //its material
                             "HDPE22");                  //its name

  HDPE_PV22 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,0*cm,25*cm),            //at (0,0,0)
                             HDPE_LV22,                      //its logical volume
                            "HDPE22",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  //The HDPE_block23

  HDPE_Box23 = new G4Box("HDPE23",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  HDPE_LV23 = new G4LogicalVolume(HDPE_Box23,                     //its shape
                              polyethylene,                      //its material
                             "HDPE23");                  //its name

  HDPE_PV23 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,0*cm,25*cm),            //at (0,0,0)
                             HDPE_LV23,                      //its logical volume
                            "HDPE23",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


 //The HDPE_block24
  

  HDPE_Box24 = new G4Box("HDPE24",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  HDPE_LV24 = new G4LogicalVolume(HDPE_Box24,                     //its shape
                              polyethylene,                      //its material
                             "HDPE24");                  //its name

  HDPE_PV24 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,10*cm,25*cm),            //at (0,0,0)
                             HDPE_LV24,                      //its logical volume
                            "HDPE24",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  //The HDPE_block25


  HDPE_Box25 = new G4Box("HDPE25",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  HDPE_LV25 = new G4LogicalVolume(HDPE_Box25,                     //its shape
                              polyethylene,                      //its material
                             "HDPE25");                  //its name

  HDPE_PV25 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,-10*cm,25*cm),            //at (0,0,0)
                             HDPE_LV25,                      //its logical volume
                            "HDPE25",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  //The lead1
  fLeadSize = 10*cm;


  Lead_Box = new G4Box("Lead",                             //its name
                   fLeadSize/2,fLeadSize/2, 5*cm/2);   //its dimensions

  Lead_LV = new G4LogicalVolume(Lead_Box,                     //its shape
                              leadMaterial,                      //its material
                             "Lead");                  //its name

  Lead_PV = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,0,3.5*cm),            //at (0,0,0)
                             Lead_LV,                      //its logical volume
                            "Lead",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


   G4VisAttributes* red = new G4VisAttributes(G4Colour::Red());

   red->SetVisibility(true);
   red->SetForceAuxEdgeVisible(true);

   Lead_LV->SetVisAttributes(red);

   //The Borated polythylene_block1 with pinhole

  BoratedSize = 30*cm;
  Borated_thickness = 3*cm;
  Borated_Box1 = new G4Box("Borated1",                             //its name
                   BoratedSize/2,  BoratedSize/2,Borated_thickness/2);   //its dimensions



  Hole = new G4Tubs("BoxHole", 0.0*cm, 0.75*cm, 1.5*cm, 0*deg, 360*deg);  // the diameter of the exit(pinhole) is 3cm. In G4 we use halfsize of the radius. 

  Hole_LV = new G4LogicalVolume(Hole,                     //its shape
                              Vacuum,                      //its material
                             "H1");                  //its name

   Borated_LV1 = new G4LogicalVolume(Borated_Box1,                     //its shape
                              b_polyethylene,                      //its material
                             "Borated1", 0,0,0);                  //its name



  Borated_PV1 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0*cm,0*cm,31.5*cm),            //at (0,0,0)
                             Borated_LV1,                      //its logical volume
                            "Borated1",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  Hole_PV = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0*cm,0*cm,0*cm),            //at (0,0,0)
                             Hole_LV,                      //its logical volume
                            "H1",                    //its name
                            Borated_LV1,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

   G4VisAttributes* green = new G4VisAttributes(G4Colour::Green());

   green->SetVisibility(true);
   green->SetForceAuxEdgeVisible(true);

   Borated_LV1->SetVisAttributes(green);


    //The lead5

  Lead_Box5 = new G4Box("Lead5",                             //its name
                   30*cm/2,30*cm/2, 3*cm/2);   //its dimensions

  Lead_LV5 = new G4LogicalVolume(Lead_Box5,                     //its shape
                              leadMaterial,                      //its material
                             "Lead5");                  //its name

  Lead_PV5 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0*cm,0*cm,34.5*cm),            //at (0,0,0)
                             Lead_LV5,                      //its logical volume
                            "Lead5",                    //its name
                            experimentalHall_log,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number
  Hole_PV2 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0*cm,0*cm,0*cm),            //at (0,0,0)
                             Hole_LV,                      //its logical volume
                            "H2",                    //its name
                             Lead_LV5,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number




   
  Lead_LV5->SetVisAttributes(red);


*/


  //Sphereball(20.00*cm);
//  ConstructValve(20.00*cm,0*cm);
  //===============PPAC DETECTOR===========================================================
  //polyethylene shield for neutron-proton conversion. =====================================
    // Create and place shield (polyethylene)
  G4double converter_x = 50.0*mm+collimatore;
  G4double converter_y = 50.0*mm+collimatore;
  G4double converter_z = 10*um;

  auto shield = new G4Box("shield", converter_x, converter_y, converter_z);
  auto lShield = new G4LogicalVolume(shield, polyethylene, "Shield");

  auto pShield = new G4PVPlacement(0,
                                      G4ThreeVector(0.0, 0.0,-1.5*mm-10*um-1.6*um+70.0*cm),
                                      lShield,
                                      "Shield",
                                      experimentalHall_log,
                                      false,
                                      0, true);
   G4VisAttributes* red = new G4VisAttributes(G4Colour::Red());

   red->SetVisibility(true);
   red->SetForceAuxEdgeVisible(true);
  
  lShield->SetVisAttributes(red);
  //PPAC
  //Particle drift volume
  G4double PPAC_drift_x = 50.0*mm+collimatore;
  G4double PPAC_drift_y = 50.0*mm+collimatore;
  G4double PPAC_drift_z = 1.5*mm;
  G4Box* PPAC_box = new G4Box("Drift Volume",PPAC_drift_x,PPAC_drift_y,PPAC_drift_z);
  PPAC_log = new G4LogicalVolume(PPAC_box,CF4,"PPAC_log",0,0,0);
  PPAC_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,70.0*cm),PPAC_log,"PPAC_Vol",experimentalHall_log,false,0,true);
  fScoringVolume_1 = PPAC_log;
 //.... using 3mm spacing between each collimator and placing the photo-sensor between each...............................................
 //Collimator Frame
  //Teflon

  G4double collf_x = collimatore/2;
  G4double collf_y = 0.5*mm;
  G4double collf_z = 1.5*mm;
  G4Box* collf_box = new G4Box("collimator",collf_x,collf_y,collf_z);
  collf_log = new G4LogicalVolume(collf_box,Teflon,"collf_log",0,0,0);
  for (int a=0;a<34;a++) 
  {
    collf_phys = new G4PVPlacement(0,G4ThreeVector(-50.0*mm-collf_x,(-49.5+a*3)*mm,0.0),collf_log,"collf_Vol",PPAC_log,false,a,true);  
  }
  for (int a=0;a<34;a++) 
  {
    collf_phys = new G4PVPlacement(0,G4ThreeVector(50.0*mm+collf_x,(-49.5+a*3)*mm,0.0),collf_log,"collf_Vol",PPAC_log,false,a+33,true);  
  }

 //Collimator Frame
  //Teflon
  G4double coll_x = 0.5*mm;
  G4double coll_y = collimatore/2;
  G4double coll_z = 1.5*mm;
  G4Box* coll_box = new G4Box("collimator",coll_x,coll_y,coll_z);
  coll_log = new G4LogicalVolume(coll_box,Teflon,"collf_log",0,0,0);
  for (int b=0;b<34;b++) 
  {
    coll_phys = new G4PVPlacement(0,G4ThreeVector((-49.5+b*3)*mm,-50.0*mm-coll_y,0.0),coll_log,"coll_Vol",PPAC_log,false,b,true);  
  }
  for (int b=0;b<34;b++) 
  {
    coll_phys = new G4PVPlacement(0,G4ThreeVector((-49.5+b*3)*mm,50.0*mm+coll_y,0.0),coll_log,"coll_Vol",PPAC_log,false,b+33,true);  
  }

  //Silicon_detector...
  //silicon
  // Define the silicon photo sensor dimensions
  G4double sensor_x  = 1.0*mm;
  G4double sensor_y  = 1.0*mm;
  G4double sensor_z = 1.5*mm; 
  G4Box* sensor_box = new G4Box("sensor", sensor_x, sensor_y, sensor_z);
  sensor_log1 = new G4LogicalVolume(sensor_box, silicon, "sensor_log1", 0, 0, 0);
  sensor_log1->SetVisAttributes(red);

  sensor_log2 = new G4LogicalVolume(sensor_box, silicon, "sensor_log2", 0, 0, 0);
  sensor_log2->SetVisAttributes(red);

  //left sensor:::::::::sensor_log1
  for (int c=0;c<33;c++) 
  {
    sensor_phys = new G4PVPlacement(0,G4ThreeVector((-50.0*mm - collimatore)+sensor_x,(-49.5+c*3+1.5)*mm,0.0),sensor_log1,"sensor_Vol1",PPAC_log,false,c,true);  
  }

  //right sensor::::::::sensor_log2

  for (int c=0;c<33;c++) 
  {
    sensor_phys = new G4PVPlacement(0,G4ThreeVector((50.0*mm + collimatore)-sensor_x,(-49.5+c*3+1.5)*mm,0.0),sensor_log2,"sensor_Vol2",PPAC_log,false,c,true);  
  }
//  fScoringVolume_1 = sensor_log;


  //sensors on the top and bottom of the PPAC........................
  // Top and Bottom Sensors
  G4double siPM_x  = 1.0*mm;
  G4double siPM_y  = 1.0*mm;
  G4Box* sensor_box2 = new G4Box("sensor2", siPM_x, siPM_y, sensor_z);
  sensor_log3 = new G4LogicalVolume(sensor_box2, silicon, "sensor_log3", 0, 0, 0);
  sensor_log3->SetVisAttributes(red);

  sensor_log4 = new G4LogicalVolume(sensor_box2, silicon, "sensor_log4", 0, 0, 0);
  sensor_log4->SetVisAttributes(red);

  // Top and Bottom Sensors
  for (int d = 0; d < 33; d++) {
      G4double sensor_x_position = (-49.5 + d * 3 + 1.5) * mm;

      // Bottom side:::::::::::::::::::sensor_log3
      sensor_phys2  = new G4PVPlacement(0, G4ThreeVector(sensor_x_position, (-50.0 * mm - collimatore) + siPM_y, 0.0),
                                    sensor_log3, "sensor_Vol3", PPAC_log, false, d , true);
  }

  for (int d = 0; d < 33; d++) {
      G4double sensor_x_position = (-49.5 + d * 3 + 1.5) * mm;

     // Top side:::::::::::::::::::::::::sensor_log4
      sensor_phys2  = new G4PVPlacement(0, G4ThreeVector(sensor_x_position, (50.0 * mm + collimatore) - siPM_y, 0.0),
                                    sensor_log4, "sensor_Vol4", PPAC_log, false, d , true);
  }
 
  //Cathode
  //Particle drift volume
  G4double cathode_x = 50.0*mm+collimatore;
  G4double cathode_y = 50.0*mm+collimatore;
  G4double cathode_z = 0.8*um;
  G4Box* cathode_box = new G4Box("cathode Volume",cathode_x,cathode_y,cathode_z);
  cathode_log = new G4LogicalVolume(cathode_box,PP,"cathode_log",0,0,0);
  cathode_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,1.5*mm+0.8*um+70.0*cm),cathode_log,"cathode_Vol",experimentalHall_log,false,0,true);

  
  //anode
  //Particle drift volume
  G4double anode_x = 50.0*mm+collimatore;
  G4double anode_y = 50.0*mm+collimatore;
  G4double anode_z = 0.8*um;
  G4Box* anode_box = new G4Box("anode Volume",anode_x,anode_y,anode_z);
  anode_log = new G4LogicalVolume(anode_box,PP,"anode_log",0,0,0);
  anode_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,-1.5*mm-0.8*um+70.0*cm),anode_log,"anode_Vol",experimentalHall_log,false,0,true);



  //Cathode-AL
  //Particle drift volume
  G4double cathodeAl_x = 50.0*mm+collimatore;
  G4double cathodeAl_y = 50.0*mm+collimatore;
  G4double cathodeAl_z = 0.05*um;
  G4Box* cathodeAl_box = new G4Box("cathodeAl Volume",cathodeAl_x,cathodeAl_y,cathodeAl_z);
  cathodeAl_log = new G4LogicalVolume(cathodeAl_box,Aluminium,"cathodeAu_log",0,0,0);
  cathodeAl_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,-0.8*um+0.05*um),cathodeAl_log,"cathodeAl_Vol",cathode_log,false,0,true);



  //Anode - AL
  //Particle drift volume
  G4double anodeAl_x = 50.0*mm+collimatore;
  G4double anodeAl_y = 50.0*mm+collimatore;
  G4double anodeAl_z = 0.05*um;
  G4Box* anodeAl_box = new G4Box("anode Volume",anodeAl_x,anodeAl_y,anodeAl_z);
  anodeAl_log = new G4LogicalVolume(anodeAl_box,Aluminium,"anodeAu_log",0,0,0);
  anodeAl_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,+0.8*um-0.05*um),anodeAl_log,"anodeAu_Vol",anode_log,false,0,true);


  // ------------------------------- Surfaces -----------------------------------------------
  // =========================================================================
       // Optical Properties
// =========================================================================

  // surface reflecting (Aluminium electrodes)
  G4OpticalSurface* oppac_Al_gas = new G4OpticalSurface("Reflecting");
  oppac_Al_gas->SetModel(unified);
  oppac_Al_gas->SetType(dielectric_metal);
  oppac_Al_gas->SetFinish(polished);
  oppac_Al_gas->SetMaterialPropertiesTable(reflectMPT);

  G4LogicalBorderSurface* oppac_1 = new G4LogicalBorderSurface("oppac_1",PPAC_phys,cathodeAl_phys,oppac_Al_gas);
  G4LogicalBorderSurface* oppac_2 = new G4LogicalBorderSurface("oppac_2",PPAC_phys,anodeAl_phys,oppac_Al_gas);
  G4LogicalBorderSurface* oppac_3 = new G4LogicalBorderSurface("oppac_3",coll_phys,cathodeAl_phys,oppac_Al_gas);
  G4LogicalBorderSurface* oppac_4 = new G4LogicalBorderSurface("oppac_4",coll_phys,anodeAl_phys,oppac_Al_gas);

 // surface assorbing  (Teflon collimators)
  G4OpticalSurface* oppac_ab = new G4OpticalSurface("Absorbing");
  oppac_ab->SetModel(unified);
  oppac_ab->SetType(dielectric_dielectric);
  oppac_ab->SetFinish(ground);
  oppac_ab->SetMaterialPropertiesTable(absorbMPT); 
  
  G4LogicalSkinSurface* oppac_10 = new G4LogicalSkinSurface("caccola_1",collf_log,oppac_ab);
  G4LogicalSkinSurface* oppac_11 = new G4LogicalSkinSurface("caccola_2",coll_log,oppac_ab);

  // CF4-Silicon detector surface

  G4LogicalSkinSurface* cf4si_1 = new G4LogicalSkinSurface("cf4si_1", sensor_log1, cf4SiSurface);
  G4LogicalSkinSurface* cf4si_2 = new G4LogicalSkinSurface("cf4si_2", sensor_log2, cf4SiSurface);
  G4LogicalSkinSurface* cf4si_3 = new G4LogicalSkinSurface("cf4si_3", sensor_log3, cf4SiSurface);
  G4LogicalSkinSurface* cf4si_4 = new G4LogicalSkinSurface("cf4si_4", sensor_log4, cf4SiSurface); 
 
}


///........................FUNCTION TO ADD THE DETECTOR TO THE SENSITIVE DETECTOR...................

void DC::ConstructSDandField()

{

    MySensitiveDetector *sensDet = new MySensitiveDetector("SensitiveDetector");

    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    SDman->AddNewDetector(sensDet);


    sensor_log1->SetSensitiveDetector(sensDet);   //left sensor
    sensor_log2->SetSensitiveDetector(sensDet);   //right sensor
    sensor_log3->SetSensitiveDetector(sensDet);   //bottom sensor
    sensor_log4->SetSensitiveDetector(sensDet);   //top sensor

}
