//-------- GEANT4
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4PVReplica.hh"
//#include "G4AssemblyVolume.hh"
//#include "G4Transform3D.hh"
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

#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4UnitsTable.hh"

#include "DC.hh"

DC::DC(double density, double lunghezza_collimatore)
 :  experimentalHall_log(0),GAS_log(0),PPAC_log(0),SiPM_log(0),collf_log(0),coll_log(0),cathode_log(0),anode_log(0),cathodeAl_log(0),anodeAl_log(0),
    experimentalHall_phys(0),GAS_phys(0),PPAC_phys(0),SiPM_phys(0),collf_phys(0),coll_phys(0),cathode_phys(0),anode_phys(0),cathodeAl_phys(0),anodeAl_phys(0)
{
  dens = density;
  collimatore = lunghezza_collimatore;
}

DC::~DC()
{;}

G4VPhysicalVolume* DC::Construct()
{
  DefineMaterials();
  ConstructLaboratory();
  SensitiveDete();
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

  //----------------------------------- vacuum ------------------------------------
  //G4Material* Vacuum = pNistManager->FindOrBuildMaterial("G4_Galactic");
  G4Material *Vacuum = new G4Material("Vacuum", 1.e-20*g/cm3, 2, kStateGas);
  Vacuum->AddElement(N, 0.755);
  Vacuum->AddElement(O, 0.245);
  //gMan->SetMaterialAppearance(vacuum,-1);

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

/*
  const G4int iNbEntries = 2;
  G4double CF4PhotonMomentum[iNbEntries] = {5*eV,4.5*eV};
  G4double CF4Scintillation_Fast[iNbEntries] = {1,0};
  G4double CF4Scintillation_Slow[iNbEntries] = {0,1};
*/

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
//  pTeflonPropertiesTable->AddProperty("SPECULARLOBECONSTANT", pdTeflonPhotonMomentum, pdTeflonSpecularLobe, TNbEntries);
//  pTeflonPropertiesTable->AddProperty("SPECULARSPIKECONSTANT", pdTeflonPhotonMomentum, pdTeflonSpecularSpike, TNbEntries);
//  pTeflonPropertiesTable->AddProperty("BACKSCATTERCONSTANT", pdTeflonPhotonMomentum, pdTeflonBackscatter, TNbEntries);
//  pTeflonPropertiesTable->AddProperty("EFFICIENCY", pdTeflonPhotonMomentum, pdTeflonEfficiency, TNbEntries);
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
//  pTeflonPropertiesTable->AddProperty("SPECULARLOBECONSTANT", pdTeflonPhotonMomentum, pdTeflonSpecularLobe, TNbEntries);
//  pTeflonPropertiesTable->AddProperty("SPECULARSPIKECONSTANT", pdTeflonPhotonMomentum, pdTeflonSpecularSpike, TNbEntries);
//  pTeflonPropertiesTable->AddProperty("BACKSCATTERCONSTANT", pdTeflonPhotonMomentum, pdTeflonBackscatter, TNbEntries);
//  pTeflonPropertiesTable->AddProperty("EFFICIENCY", pdTeflonPhotonMomentum, pdTeflonEfficiency, TNbEntries);

}

void DC::ConstructLaboratory()
{
  G4Material *Vacuum = G4Material::GetMaterial("Vacuum");
  G4Material *CF4 = G4Material::GetMaterial("CF4");
  G4Material *Teflon = G4Material::GetMaterial("Teflon");
  G4Material *PP = G4Material::GetMaterial("G4_POLYPROPYLENE");
  //G4Material *Aluminium = G4Material::GetMaterial("Aluminium");

  //	========================================================================================================================
  //======================================================================================================================
  // Experimental hall (world volume) --------------
  G4double expHall_x = 5.0*m;
  G4double expHall_y = 5.0*m;
  G4double expHall_z = 5.0*m;
  G4Box* experimentalHall_box = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);
  experimentalHall_log = new G4LogicalVolume(experimentalHall_box,Vacuum,"expHall_log");
  experimentalHall_phys = new G4PVPlacement(0,G4ThreeVector(),experimentalHall_log,"expHall",0,false,0);
  experimentalHall_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  //polyethylene shield for neutron-proton conversion. =====================================
    // Create and place shield (polyethylene)

  G4double converter_x = 50.0*mm+collimatore;
  G4double converter_y = 50.0*mm+collimatore;
  G4double converter_z = 10*um;

  auto shield = new G4Box("shield", converter_x, converter_y, converter_z);
  auto lShield = new G4LogicalVolume(shield, polyethylene, "Shield");

  auto pShield = new G4PVPlacement(0,
                                      G4ThreeVector(0.0, 0.0,-1.5*mm-10*um-1.6*um+112*cm),
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
  PPAC_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,112.0*cm),PPAC_log,"PPAC_Vol",experimentalHall_log,false,0,true);
  fScoringVolume_1 = PPAC_log;
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

  //Cathode
  //Particle drift volume
  G4double cathode_x = 50.0*mm+collimatore;
  G4double cathode_y = 50.0*mm+collimatore;
  G4double cathode_z = 0.8*um;
  G4Box* cathode_box = new G4Box("cathode Volume",cathode_x,cathode_y,cathode_z);
  cathode_log = new G4LogicalVolume(cathode_box,PP,"cathode_log",0,0,0);
  cathode_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,1.5*mm+0.8*um+112.0*cm),cathode_log,"cathode_Vol",experimentalHall_log,false,0,true);

  
  //anode
  //Particle drift volume
  G4double anode_x = 50.0*mm+collimatore;
  G4double anode_y = 50.0*mm+collimatore;
  G4double anode_z = 0.8*um;
  G4Box* anode_box = new G4Box("anode Volume",anode_x,anode_y,anode_z);
  anode_log = new G4LogicalVolume(anode_box,PP,"anode_log",0,0,0);
  anode_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,-1.5*mm-0.8*um+112.0*cm),anode_log,"anode_Vol",experimentalHall_log,false,0,true);



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
  // ------------------------- Optical Propertiers ------------------------------------------
 
  // surface reflecting
  G4OpticalSurface* oppac_Al_gas = new G4OpticalSurface("Reflecting");
  oppac_Al_gas->SetModel(unified);
  oppac_Al_gas->SetType(dielectric_metal);
  oppac_Al_gas->SetFinish(polished);
  G4LogicalBorderSurface* oppac_1 = new G4LogicalBorderSurface("oppac_1",PPAC_phys,cathodeAl_phys,oppac_Al_gas);
  G4LogicalBorderSurface* oppac_2 = new G4LogicalBorderSurface("oppac_2",PPAC_phys,anodeAl_phys,oppac_Al_gas);
  G4LogicalBorderSurface* oppac_3 = new G4LogicalBorderSurface("oppac_3",coll_phys,cathodeAl_phys,oppac_Al_gas);
  G4LogicalBorderSurface* oppac_4 = new G4LogicalBorderSurface("oppac_4",coll_phys,anodeAl_phys,oppac_Al_gas);

 // surface assorbing
  G4OpticalSurface* oppac_ab = new G4OpticalSurface("Absorbing");
  oppac_ab->SetModel(unified);
  oppac_ab->SetType(dielectric_dielectric);
//  oppac_ab->SetFinish(polished);
//  oppac_ab->SetFinish(polishedbackpainted);
  oppac_ab->SetFinish(ground);
  G4LogicalSkinSurface* oppac_10 = new G4LogicalSkinSurface("caccola_1",collf_log,oppac_ab);
  G4LogicalSkinSurface* oppac_11 = new G4LogicalSkinSurface("caccola_2",coll_log,oppac_ab);


 
}


void DC::SensitiveDete()
{
/*
  // ------------------------------------------------------------------------------------------
  // sensitive detectors -----------------------------------------------------
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  SDman->SetVerboseLevel(0);
  //--------------------------------------------------------------------------------------------
  //Define Multi-Detector and Register 
  //--------------------------------------------------------------------------------------------
  G4MultiFunctionalDetector* det = new G4MultiFunctionalDetector("IonPro");
  SDman->AddNewDetector(det);
  PPAC_log->SetSensitiveDetector(det);
  //--------------------------------------------------------------------------------------------
  //Primitive Score
  //--------------------------------------------------------------------------------------------
  G4VPrimitiveScorer* primitive;
  primitive = new G4PSEnergyDeposit("nproton");
  //--------------------------------------------------------------------------------------------
  //Filter
  //--------------------------------------------------------------------------------------------
  G4String filterName,particleName;
  //-- neutron filter
  G4SDParticleFilter* epFilter = new G4SDParticleFilter(filterName="epFilter");
  epFilter->add(particleName="alpha");
  //epFilter->add(particleName="neutron");
  primitive->SetFilter(epFilter);
  //epFilter->show();
  //-- Register Sensitive Detector Run Action 
  det->RegisterPrimitive(primitive);*/
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
}

