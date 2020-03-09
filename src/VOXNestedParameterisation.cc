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


#include "../include/VOXNestedParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"

#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"
#include "G4Material.hh"


VOXNestedParameterisation::VOXNestedParameterisation(VOXModelImport* _voxelPhantom)
:voxelPhantom(_voxelPhantom)
{
	blankAtt = new G4VisAttributes(G4VisAttributes::Invisible);
	blankAtt->SetVisibility(FALSE);

	fNx = voxelPhantom->GetVoxelResolution(0);
	fNy = voxelPhantom->GetVoxelResolution(1);
	fNz = voxelPhantom->GetVoxelResolution(2);

	fVoxelHalfLengthX = voxelPhantom->GetVoxelSize(0) * 0.5;
	fVoxelHalfLengthY = voxelPhantom->GetVoxelSize(1) * 0.5;
	fVoxelHalfLengthZ = voxelPhantom->GetVoxelSize(2) * 0.5;
}

VOXNestedParameterisation::~VOXNestedParameterisation()
{
}


G4Material* VOXNestedParameterisation::
ComputeMaterial(G4VPhysicalVolume* physVol, const G4int iz,
                const G4VTouchable* parentTouch)
{
    // protection for initialization and vis at idle state
    //
	if(parentTouch==0) return (voxelPhantom->GetVoxelMaterial(0));

	G4int ix = parentTouch->GetReplicaNumber(0);
	G4int iy = parentTouch->GetReplicaNumber(1);

	G4int voxelData = voxelPhantom->GetVoxelData(ix,iy,iz);

    return voxelPhantom->GetVoxelMaterial(voxelData);


}

unsigned int VOXNestedParameterisation::
GetMaterialIndex( unsigned int copyNo ) const
{
    return copyNo;
}


G4int VOXNestedParameterisation::GetNumberOfMaterials() const
{
    return voxelPhantom->GetVoxelMaterialSize();
}

G4Material* VOXNestedParameterisation::GetMaterial(G4int i) const
{
	if(i >= voxelPhantom->GetVoxelMaterialSize()) i=0;
	G4int idx = voxelPhantom->GetMaterialIndex(i);
	return voxelPhantom->GetVoxelMaterial(idx);
}


//
// Transformation of voxels.
//
void VOXNestedParameterisation::
ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
    // Position of voxels.
    // x and y positions are already defined in DetectorConstruction by using
    // replicated volume. Here only we need to define is z positions of voxels.
    physVol->SetTranslation(G4ThreeVector(0.,0.,(2.*static_cast<double>(copyNo)+1.)*fVoxelHalfLengthZ - fVoxelHalfLengthZ*fNz));
}

//
// Dimensions are always same in this RE02 example.
//
void VOXNestedParameterisation::
ComputeDimensions( G4Box& box, const G4int, const G4VPhysicalVolume* ) const
{
    box.SetXHalfLength(fVoxelHalfLengthX);
    box.SetYHalfLength(fVoxelHalfLengthY);
    box.SetZHalfLength(fVoxelHalfLengthZ);
}
