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
// TETParameterisation.cc
// \file   MRCP_GEANT4/External/src/TETParameterisation.cc
// \author Haegin Han
//

#include "TETParameterisation.hh"
#include "G4LogicalVolume.hh"
#include "G4VisExecutive.hh"
#include "G4RunManager.hh"

TETParameterisation::TETParameterisation(TETModelImport* _tetData)
: G4VPVParameterisation(), tetData(_tetData)
{
	// initialise visAttMap which contains G4VisAttributes* for each organ
	auto colourMap =  tetData->GetColourMap();
	for(auto colour : colourMap){
		visAttMap[colour.first] = new G4VisAttributes(colour.second);
	}

	if(colourMap.size()) isforVis = true;
	else                 isforVis = false;
}

TETParameterisation::~TETParameterisation()
{}

G4VSolid* TETParameterisation::ComputeSolid(
    		       const G4int copyNo, G4VPhysicalVolume* )
{
	// return G4Tet*
	return tetData->GetTetrahedron(copyNo);
}

void TETParameterisation::ComputeTransformation(
                   const G4int,G4VPhysicalVolume*) const
{}

G4Material* TETParameterisation::ComputeMaterial(const G4int copyNo,
                                                 G4VPhysicalVolume* phy,
                                                 const G4VTouchable* )
{
   // set the colour for each organ if visualization is required
	if(isforVis){
		G4int idx = tetData->GetMaterialIndex(copyNo);
		if(visAttMap.find(idx)==visAttMap.end()) idx = visAttMap.begin()->first;
		phy->GetLogicalVolume()->SetVisAttributes(visAttMap[idx]);
		phy->GetLogicalVolume()->SetMaterial(tetData->GetMaterial(tetData->GetMaterialIndex(copyNo)));
	}

	// return the material data for each material index
	return tetData->GetMaterial(tetData->GetMaterialIndex(copyNo));
}


