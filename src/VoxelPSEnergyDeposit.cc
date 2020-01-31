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
// TETPSEnergyDeposit.cc
// \file   MRCP_GEANT4/External/src/TETPSEnergyDeposit.cc
// \author Haegin Han
//

#include "VoxelPSEnergyDeposit.hh"

VoxelPSEnergyDeposit::VoxelPSEnergyDeposit(G4String name, VOXModelImport* _voxData)
  :G4PSEnergyDeposit(name), voxData(_voxData)
{}

VoxelPSEnergyDeposit::~VoxelPSEnergyDeposit()
{}

G4int VoxelPSEnergyDeposit::GetIndex(G4Step* aStep)
{
	const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
	G4int ix = touchable->GetReplicaNumber(1);
	G4int iy = touchable->GetReplicaNumber(2);
	G4int iz = touchable->GetReplicaNumber(0);

	G4int organID = voxData->GetVoxelData(ix,iy,iz);

	return organID;
}
