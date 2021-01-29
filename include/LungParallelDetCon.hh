//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// The code was written by : hhg
//
// Department of Nuclear Engineering, Hanyang University
// 17 Haengdang, Seongdong, Seoul 133-791, Korea
// Tel: +82-2-2220-4053
// Fax: +82-2-2220-4059
//
#ifndef LungParallelWorldConstruction_h
#define LungParallelWorldConstruction_h 1

#include "G4VUserParallelWorld.hh"
#include "globals.hh"
#include "TETModelImport.hh"
#include "LungDetMessenger.hh"

#include <fstream>
#include <vector>
#include <utility>
#include <map>

typedef std::pair<G4ThreeVector,G4bool> BB;

class G4LogicalVolume;
class G4VPhysicalVolume;

class LungParallelDetCon : public G4VUserParallelWorld
{
  public:
    LungParallelDetCon(G4String parallelWorldName, TETModelImport* tetData);
    virtual ~LungParallelDetCon();

  public:
    virtual void Construct();
    virtual void ConstructSD();
    std::pair<G4ThreeVector,G4ThreeVector> GetBbox() {return bbox;}

  private:
    void ImportInfoData();
    void SetBBdiameters();
    void PrintLungData();
    void ConstructBBUnit(std::vector<BB> fBB_vec, G4int mid_gen);
    void SetUpDetectors();

    G4bool fConstructed;
    TETModelImport* tetData;
    G4double bb_diam[17][10];
    G4ThreeVector box_center;
    G4LogicalVolume* bBoxLogical;
    static G4ThreadLocal G4bool fSDConstructed;
    std::vector<BB> bb_vec;
    std::vector<G4VPhysicalVolume*> phys_vec;
    std::map<int, G4Material*> mat_map;
    std::pair<G4ThreeVector,G4ThreeVector> bbox;

    //for messenger
public:
    void SetVolChkName(G4String fileN) {volChkName=fileN;}
    void SetBBbasVol(G4double vol)     {tetData->SetBB_basal_vol(vol);}
    void SetBBsecVol(G4double vol)     {tetData->SetBB_secretory_vol(vol);}
    void SetbbsecVol(G4double vol)     {tetData->Setbb_secretory_vol(vol);}

private:
    LungDetMessenger* fMessenger;
    G4String volChkName;
};

#endif


