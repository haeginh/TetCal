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
/// \file RE06/include/SeedParallel.hh
/// \brief Definition of the SeedParallel class
//
// 

#ifndef SeedParallel_h
#define SeedParallel_h 1

#include "G4VUserParallelWorld.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include <vector>

class G4Tet;

class SeedParallel : public G4VUserParallelWorld
{
  public:
    SeedParallel(G4String worldName);
    virtual ~SeedParallel();

    virtual void Construct();
    virtual void ConstructSD();
    G4int GetNumOfTet() {return tetVec.size();}
    G4Tet* GetTet(G4int i) {return tetVec[i];}
    G4int GetId(G4int i) {return matVec[i];}
    G4ThreeVector GetCenter() {return (bBox_max+bBox_min)*0.5;}

  private:
    G4bool fConstructed;
    void SetupGeometry();
    void ReadTet();
    std::vector<G4Tet*> tetVec; 
    std::vector<G4int> matVec;
    std::map<G4int, G4Material*> materials;
    G4ThreeVector bBox_max, bBox_min;
};


#endif

