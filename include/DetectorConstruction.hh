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
//	*Yeon Soo Yeom, yeonsoo.yeom@nih.gov
//	*Choonsik Lee, choonsik.lee@nih.gov
//
// Radiation Epidemiology Branch, DCEG/NCI/NHI
// 9609 Medical Center Dr, Rociville, MD 20850
// Tel: +1-240-276-532
// ********************************************************************

#ifndef DETECTORCONSTRUCTION_HH
#define DETECTORCONSTRUCTION_HH

#include <fstream>
#include <sstream>

#include "globals.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VisAttributes.hh"

#include "VoxelNestedParameterisation.hh"
#include "DRFScorer.hh"
#include "ImportVoxelPhantom.hh"

class DetectorConstruction : public G4VUserDetectorConstruction {

public:

	DetectorConstruction(ImportVoxelPhantom*, G4String);
	virtual ~DetectorConstruction();

	virtual G4VPhysicalVolume* Construct();
	virtual void ConstructSDandField();
	
	G4Transform3D ReadColliAndDAPInfo();
	void ResetColli(G4String fileName);
  
private:
	G4VisAttributes* blankAtt;

	ImportVoxelPhantom* voxelPhantom;

	G4Material* air;

	// world & container logical and physical volumes
	G4VPhysicalVolume* fpWorldPhysical;
	G4LogicalVolume* logicVoxel;
	// G4LogicalVolume* logDAP;
	
	// Data members
	G4ThreeVector fphantomSize;   // Size of voxel phantom
	G4ThreeVector fvoxelSize;     // voxel size
	G4int         fNx,fNy,fNz;    // Number of segmentation of voxel phantom
	G4double      fVoxelHalfLengthX, fVoxelHalfLengthY, fVoxelHalfLengthZ;    // Number of segmentation of voxel phantom
	
	G4double BeamPortHalfSizeX, BeamPortHalfSizeY, BeamPortHalfSizeZ;
	G4double CollimatorCenterX, CollimatorCenterY, CollimatorCenterZ;
	G4double tableDen;
	G4ThreeVector tableHalfSize, tablePos;
	// G4double DAPHalfSizeX, DAPHalfSizeY, DAPHalfSizeZ;
	// G4double DAPCenterX, DAPCenterY, DAPCenterZ;
	G4double DirectionX, DirectionY, DirectionZ;
	G4double RotAngle;
	G4ThreeVector RotAxis;
	G4RotationMatrix* Rot;
	G4Transform3D transform;
	
	G4String inp_name;

	G4VSolid* solCollimator;
	G4VSolid* solBeamPort;
	G4VSolid* solCollimator_F;
	G4LogicalVolume* logCollimator;
	G4VPhysicalVolume* phyCollimator;

	
};

#endif

