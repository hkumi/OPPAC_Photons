#ifndef __DC_H__
#define __DC_H__

#include "G4Material.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

class DC : public G4VUserDetectorConstruction
{
  public:
    DC(double density, double lunghezza_collimatore);
    virtual ~DC();

  public:
    virtual G4VPhysicalVolume* Construct();
    G4LogicalVolume *GetScoringVolume() const {return fScoringVolume_1;}
  private:
    // Logical volumes
    G4LogicalVolume* experimentalHall_log;
    G4LogicalVolume* GAS_log;
    G4LogicalVolume* PPAC_log;
    G4LogicalVolume* SiPM_log;
    G4LogicalVolume* collf_log;
    G4LogicalVolume* coll_log;
    G4LogicalVolume* cathode_log;
    G4LogicalVolume* anode_log;
    G4LogicalVolume* cathodeAl_log;
    G4LogicalVolume* anodeAl_log;
    // Physical volumes
    G4VPhysicalVolume* experimentalHall_phys;
    G4VPhysicalVolume* GAS_phys;
    G4VPhysicalVolume* PPAC_phys;
    G4VPhysicalVolume* SiPM_phys;
    G4VPhysicalVolume* collf_phys;
    G4VPhysicalVolume* coll_phys;
    G4VPhysicalVolume* cathode_phys;
    G4VPhysicalVolume* anode_phys;
    G4VPhysicalVolume* cathodeAl_phys;
    G4VPhysicalVolume* anodeAl_phys;

   //Scorer.
    G4LogicalVolume   *fScoringVolume_1;

 private:
    void DefineMaterials();
    void ConstructLaboratory();
    void SensitiveDete();
    double dens;
    double collimatore;
    

 private:
    G4double LayerThickness;
    G4int NbOfLayers;
    G4Material  *b_polyethylene,  *polyethylene, *NaI,*mat_graphite;
    G4Material  *leadMaterial,*Aluminium,*PP;
    G4Element  *Na, *I, *C,*Al;
    
  

};

#endif
