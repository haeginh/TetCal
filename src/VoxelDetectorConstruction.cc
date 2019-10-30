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
// The code was written by :
//	*Jong Hwi Jeong  jonghwi@hanyang.ac.kr
//      *Chan Hyeong Kim chkim@hanyang.ac.kr
//
// Department of Nuclear Engineering, Hanyang University
// 17 Haengdang, Seongdong, Seoul 133-791, Korea
// Tel: +82-2-2220-4057
// Fax: +82-2-2220-4059
//
// ********************************************************************
#include "VoxelDetectorConstruction.hh"

#include "globals.hh"

#include "G4PhysicalConstants.hh"
#include "G4NistManager.hh"
#include "G4LogicalVolume.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"

VoxelDetectorConstruction::VoxelDetectorConstruction(ImportVoxelPhantom* _voxelPhantom)
:fpWorldPhysical(0), logicVoxel(0), voxelPhantom(_voxelPhantom)
{

	blankAtt = new G4VisAttributes;
	blankAtt->SetVisibility(FALSE);

	fNx = voxelPhantom->GetVoxelResolution(0);
	fNy = voxelPhantom->GetVoxelResolution(1);
	fNz = voxelPhantom->GetVoxelResolution(2);

	fVoxelHalfLengthX = voxelPhantom->GetVoxelSize(0) * 0.5;
	fVoxelHalfLengthY = voxelPhantom->GetVoxelSize(1) * 0.5;
	fVoxelHalfLengthZ = voxelPhantom->GetVoxelSize(2) * 0.5;

	G4cout << fNx << "\t"<< fNy << "\t"<< fNz << "\t"<<G4endl;
	G4cout << fVoxelHalfLengthX << "\t"<< fVoxelHalfLengthY << "\t"<< fVoxelHalfLengthZ << "\t"<<G4endl;


}

VoxelDetectorConstruction::~VoxelDetectorConstruction()
{}

G4VPhysicalVolume* VoxelDetectorConstruction::Construct()
{


	G4Material* VACUUM = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
	// World
	G4VSolid* sol_World = new G4Box("World", 10*m, 10*m, 10*m);
	//G4Material* AIR = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
	G4LogicalVolume* lv_World = new G4LogicalVolume(sol_World, VACUUM, "World");
	fpWorldPhysical =
		new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), lv_World, "World", 0, false, 1);

	// Phantom Geometry
	// 1st - Container for Phantom
	G4String ContainerName("Container");
	G4VSolid* solContainer = new G4Box(ContainerName,
						fVoxelHalfLengthX * fNx,
						fVoxelHalfLengthY * fNy,
						fVoxelHalfLengthZ * fNz);

	G4LogicalVolume* logContainer = new G4LogicalVolume(solContainer,VACUUM,ContainerName);
	new G4PVPlacement(0, G4ThreeVector(0,0,0), logContainer, ContainerName, lv_World, false, 1);


	// Replication of Water Phantom Volume.
	// 2nd - Y Slice
	G4String yRepName("RepY");
	G4VSolid* solYRep = new G4Box(yRepName,fNx*fVoxelHalfLengthX, fVoxelHalfLengthY, fNz*fVoxelHalfLengthZ);
	G4LogicalVolume* logYRep = new G4LogicalVolume(solYRep,VACUUM,yRepName);
	new G4PVReplica(yRepName,logYRep,logContainer,kYAxis,fNy,fVoxelHalfLengthY*2.);

	logYRep->SetVisAttributes(blankAtt);

	// 3rd - X Slice
	G4String xRepName("RepX");
	G4VSolid* solXRep = new G4Box(xRepName,fVoxelHalfLengthX, fVoxelHalfLengthY, fNz*fVoxelHalfLengthZ);
	G4LogicalVolume* logXRep = new G4LogicalVolume(solXRep,VACUUM,xRepName);
	new G4PVReplica(xRepName,logXRep,logYRep,kXAxis,fNx,fVoxelHalfLengthX*2.);


	logXRep->SetVisAttributes(blankAtt);

	// Voxel solid and logical volumes
	// 4th - Z Slice
	G4VSolid* solVoxel = new G4Box("phantom",fVoxelHalfLengthX, fVoxelHalfLengthY,fVoxelHalfLengthZ);
	logicVoxel = new G4LogicalVolume(solVoxel,VACUUM,"VoxelDetector");
	logicVoxel->SetVisAttributes(blankAtt);

	//
	// Parameterisation for transformation of voxels.
	//  (voxel size is fixed in this example.
	//    e.g. nested parameterisation handles material
	//    and transfomation of voxels.)

	VoxelNestedParameterisation* param =
			new VoxelNestedParameterisation(voxelPhantom);

	new G4PVParameterised("phantom",    // their name
						  logicVoxel, // their logical volume
						  logXRep,      // Mother logical volume
						  kZAxis,       // Are placed along this axis
						  //kUndefined,	  // Are placed along this axis
						  fNz,      // Number of cells
						  param);       // Parameterisation.


	// Visualization
/*	G4VisAttributes* va_World = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	va_World->SetForceWireframe(true);
	lv_World->SetVisAttributes(va_World);
*/
	return fpWorldPhysical;
}

void VoxelDetectorConstruction::ConstructSDandField()
{
	G4SDManager* pSDman = G4SDManager::GetSDMpointer();
	G4String phantomSDname = "PhantomSD";

	G4MultiFunctionalDetector* MFDet = new G4MultiFunctionalDetector(phantomSDname);
	pSDman->AddNewDetector( MFDet );

	G4VPrimitiveScorer* scorer1 = new VoxelPSEnergyDeposit("eDep", voxelPhantom);

	MFDet->RegisterPrimitive(scorer1);

	SetSensitiveDetector(logicVoxel, MFDet);
}















