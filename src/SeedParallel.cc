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
/// \file RE06/src/SeedParallel.cc
/// \brief Implementation of the SeedParallel class
//
// 

#include "SeedParallel.hh"

#include "G4Tet.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "TETParameterisation_EyePlaque.hh"

#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"

#include <fstream>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SeedParallel::SeedParallel(G4String worldName)
:G4VUserParallelWorld(worldName),
 fConstructed(false), bBox_max(-__DBL_MAX__, -__DBL_MAX__, -__DBL_MAX__), bBox_min(__DBL_MAX__, __DBL_MAX__, __DBL_MAX__)
{
  ReadTet();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SeedParallel::~SeedParallel()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SeedParallel::Construct()
{
  if(!fConstructed)
  { 
    fConstructed = true;
    SetupGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SeedParallel::ConstructSD()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SeedParallel::SetupGeometry()
{
  //material
  std::map<G4int, G4Material*> materials;
  materials[1]=G4NistManager::Instance()->FindOrBuildMaterial("G4_Ti");
  materials[2]=G4NistManager::Instance()->FindOrBuildMaterial("G4_Ar");
  materials[3]=G4NistManager::Instance()->FindOrBuildMaterial("G4_Ag");
  materials[4]=new G4Material("coating", 6.003*g/cm3, 2, kStateSolid);
  materials[4]->AddElement(G4NistManager::Instance()->FindOrBuildElement(47), 0.5);
  materials[4]->AddElement(G4NistManager::Instance()->FindOrBuildElement(53), 0.5);
  materials[5]=G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
  
  //     
  // World
  //
  G4VPhysicalVolume* ghostWorld = GetWorld();
  G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();
  worldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
  
  G4ThreeVector halfSize = (bBox_max-bBox_min)*0.5;
  G4ThreeVector center = (bBox_max+bBox_min)*0.5;
	G4Box* containerSolid = new G4Box("phantomBox", halfSize.x() + 0.1*cm,
										                halfSize.y() + 0.1*cm,
										                halfSize.z() + 0.1*cm);

	auto containerLogic = new G4LogicalVolume(containerSolid, 0, "eyePlaque");
  new G4PVPlacement(0, center, containerLogic, "eyePlaque",
			              worldLogical, false, 0);
  G4VSolid* tetraSolid = new G4Tet("TetSolid",
			                    G4ThreeVector(),
			                    G4ThreeVector(1.*cm,0,0),
			                    G4ThreeVector(0,1.*cm,0),
			                    G4ThreeVector(0,0,1.*cm));
  auto tetLogic = new G4LogicalVolume(tetraSolid, 0, "TetLogic");

	// physical volume (phantom) constructed as parameterised geometry
	new G4PVParameterised("eyePlaque",tetLogic,containerLogic,
			                  kUndefined,tetVec.size(),
						            new TETParameterisation_EyePlaque(tetVec, matVec, materials));
}

void SeedParallel::ReadTet(){
  std::ifstream ifs("../EyePlaque1.node");
  G4int num, tmp;
  ifs>>num>>tmp>>tmp>>tmp;
  G4ThreeVector node;
  std::vector<G4ThreeVector> nodeVec;
  for(G4int i=0;i<num;i++)
  {
    ifs>>tmp>>node;
    node *= cm;
    if(node.getX()>bBox_max.getX()) bBox_max.setX(node.getX());
    if(node.getY()>bBox_max.getY()) bBox_max.setY(node.getY());
    if(node.getZ()>bBox_max.getZ()) bBox_max.setZ(node.getZ());
    if(node.getX()<bBox_min.getX()) bBox_min.setX(node.getX());
    if(node.getY()<bBox_min.getY()) bBox_min.setY(node.getY());
    if(node.getZ()<bBox_min.getZ()) bBox_min.setZ(node.getZ());
    nodeVec.push_back(node);
  }
  ifs.close();

  std::ifstream ifs2("../EyePlaque2.node");
  ifs2>>num>>tmp>>tmp>>tmp;
  std::vector<G4ThreeVector> nodeVec2;
  for(G4int i=0;i<num;i++)
  {
    ifs2>>tmp>>node;
    node *= cm;
    if(node.getX()>bBox_max.getX()) bBox_max.setX(node.getX());
    if(node.getY()>bBox_max.getY()) bBox_max.setY(node.getY());
    if(node.getZ()>bBox_max.getZ()) bBox_max.setZ(node.getZ());
    if(node.getX()<bBox_min.getX()) bBox_min.setX(node.getX());
    if(node.getY()<bBox_min.getY()) bBox_min.setY(node.getY());
    if(node.getZ()<bBox_min.getZ()) bBox_min.setZ(node.getZ());
    nodeVec2.push_back(node);
  }
  ifs2.close();

  G4ThreeVector center = (bBox_max+bBox_min)*0.5;
  std::ifstream ifsEle("../EyePlaque1.ele");
  G4int id, a, b, c, d;
  ifsEle>>num>>tmp>>tmp;
  for(G4int i=0;i<num;i++){
    ifsEle>>tmp>>a>>b>>c>>d>>id;
    tetVec.push_back(new G4Tet("tet", nodeVec[a]-center, nodeVec[b]-center, nodeVec[c]-center, nodeVec[d]-center));
    matVec.push_back(id);
  }
  ifsEle.close();
  std::ifstream ifsEle2("../EyePlaque2.ele");
  ifsEle2>>num>>tmp>>tmp;
  for(G4int i=0;i<num;i++){
    ifsEle2>>tmp>>a>>b>>c>>d>>id;
    tetVec.push_back(new G4Tet("tet", nodeVec2[a]-center, nodeVec2[b]-center, nodeVec2[c]-center, nodeVec2[d]-center));
    matVec.push_back(id);
  }
  ifsEle2.close();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
