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
// TETDetectorConstruction.cc
// \file   MRCP_GEANT4/External/src/TETDetectorConstruction.cc
// \author Haegin Han
//

#include "TETDetectorConstruction.hh"
#include "TETPSTrackLength.hh"
#include "G4VisAttributes.hh"

TETDetectorConstruction::TETDetectorConstruction(TETModelImport* _tetData)
:worldPhysical(0), container_logic(0), tetData(_tetData), tetLogic(0)
{
	// initialisation of the variables for phantom information
	phantomSize     = tetData -> GetPhantomSize();
	phantomBoxMin   = tetData -> GetPhantomBoxMin();
	phantomBoxMax   = tetData -> GetPhantomBoxMax();
	nOfTetrahedrons = tetData -> GetNumTetrahedron();
}

TETDetectorConstruction::~TETDetectorConstruction()
{
	delete tetData;
}

G4VPhysicalVolume* TETDetectorConstruction::Construct()
{
	SetupWorldGeometry();
	ConstructPhantom();
	PrintPhantomInformation();
	return worldPhysical;
}

void TETDetectorConstruction::SetupWorldGeometry()
{
	// Define the world box (size: 10*10*10 m3)
	//
	G4double worldXYZ = 10. * m;
	G4Material* vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

	G4VSolid* worldSolid
	  = new G4Box("worldSolid", worldXYZ/2, worldXYZ/2, worldXYZ/2);

	G4LogicalVolume* worldLogical
	  = new G4LogicalVolume(worldSolid,vacuum,"worldLogical");

	worldPhysical
	  = new G4PVPlacement(0,G4ThreeVector(), worldLogical,"worldPhysical", 0, false,0,false);

	// Define the phantom container (10-cm margins from the bounding box of phantom)
	//
	G4Box* containerSolid = new G4Box("phantomBox", phantomSize.x()/2 + 10.*cm,
										            phantomSize.y()/2 + 10.*cm,
										            phantomSize.z()/2 + 10.*cm);

	container_logic = new G4LogicalVolume(containerSolid, vacuum, "phantomLogical");

	new G4PVPlacement(0, G4ThreeVector(), container_logic, "PhantomPhysical",
			          worldLogical, false, 0);
	container_logic->SetOptimisation(TRUE);
	container_logic->SetSmartless( 0.5 ); // for optimization (default=2)
}

void TETDetectorConstruction::ConstructPhantom()
{
	// Define the tetrahedral mesh phantom as a parameterised geometry
	//
	// solid and logical volume to be used for parameterised geometry
	G4VSolid* tetraSolid = new G4Tet("TetSolid",
			                    G4ThreeVector(),
			                    G4ThreeVector(1.*cm,0,0),
			                    G4ThreeVector(0,1.*cm,0),
			                    G4ThreeVector(0,0,1.*cm));

	G4Material* vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
	tetLogic = new G4LogicalVolume(tetraSolid, vacuum, "TetLogic");

	// physical volume (phantom) constructed as parameterised geometry
	new G4PVParameterised("wholePhantom",tetLogic,container_logic,
			              kUndefined,tetData->GetNumTetrahedron(),
						  new TETParameterisation(tetData));
}

void TETDetectorConstruction::ConstructSDandField()
{
	// Define detector (Phantom SD) and scorer (eDep)
	//
	G4SDManager* pSDman = G4SDManager::GetSDMpointer();
	G4String phantomSDname = "PhantomSD";

	// MultiFunctional detector
	G4MultiFunctionalDetector* MFDet = new G4MultiFunctionalDetector(phantomSDname);
	pSDman->AddNewDetector( MFDet );

	// scorer for energy depositon in each organ
	MFDet->RegisterPrimitive(new TETPSTrackLength("tLengh"));

	// attach the detector to logical volume for parameterised geometry (phantom geometry)
	SetSensitiveDetector(tetLogic, MFDet);
}

void TETDetectorConstruction::PrintPhantomInformation()
{
	// print brief information on the imported phantom
	G4cout<< G4endl;
	G4cout.precision(3);
	G4cout<<"   Phantom name               "<<tetData->GetPhantomName() << " TET phantom"<<G4endl;
	G4cout<<"   Phantom size               "<<phantomSize.x()<<" * "<<phantomSize.y()<<" * "<<phantomSize.z()<<" mm3"<<G4endl;
	G4cout<<"   Phantom box position (min) "<<phantomBoxMin.x()<<" mm, "<<phantomBoxMin.y()<<" mm, "<<phantomBoxMin.z()<<" mm"<<G4endl;
	G4cout<<"   Phantom box position (max) "<<phantomBoxMax.x()<<" mm, "<<phantomBoxMax.y()<<" mm, "<<phantomBoxMax.z()<<" mm"<<G4endl;
	G4cout<<"   Number of tetrahedrons     "<<nOfTetrahedrons<<G4endl<<G4endl;
}
