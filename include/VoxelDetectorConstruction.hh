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

#ifndef DETECTORCONSTRUCTION_HH
#define DETECTORCONSTRUCTION_HH

#include <fstream>
#include <sstream>

#include "globals.hh"
#include "G4Box.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VisAttributes.hh"

#include "VoxelNestedParameterisation.hh"
#include "VoxelPSEnergyDeposit.hh"
#include "ImportVoxelPhantom.hh"

class VoxelDetectorConstruction : public G4VUserDetectorConstruction {

public:

	VoxelDetectorConstruction(VOXModelImport* );
	virtual ~VoxelDetectorConstruction();

	virtual G4VPhysicalVolume* Construct();
	virtual void ConstructSDandField();

private:
	G4VisAttributes* blankAtt;

	VOXModelImport* voxelPhantom;

	// world & container logical and physical volumes
	G4VPhysicalVolume* fpWorldPhysical;
	G4LogicalVolume* logicVoxel;
	// Data members
	G4ThreeVector fphantomSize;   // Size of voxel phantom
	G4ThreeVector fvoxelSize;     // voxel size
	G4int         fNx,fNy,fNz;    // Number of segmentation of voxel phantom
	G4double      fVoxelHalfLengthX, fVoxelHalfLengthY, fVoxelHalfLengthZ;    // Number of segmentation of voxel phantom
};

#endif

