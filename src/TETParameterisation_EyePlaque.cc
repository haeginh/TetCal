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
// TETParameterisation_EyePlaque.cc
// \file   MRCP_GEANT4/External/src/TETParameterisation_EyePlaque.cc
// \author Haegin Han
//

#include "TETParameterisation_EyePlaque.hh"
#include "G4LogicalVolume.hh"
#include "G4VisExecutive.hh"
#include "G4RunManager.hh"

TETParameterisation_EyePlaque::TETParameterisation_EyePlaque(std::vector<G4Tet*> _tetVec, std::vector<G4int> _idVec, std::map<G4int, G4Material*> _matMap)
: G4VPVParameterisation(), tetVec(_tetVec), idVec(_idVec), matMap(_matMap)
{
	// initialise visAttMap which contains G4VisAttributes* for each organ
	visAttMap[1] = new G4VisAttributes(G4Color(1., 0., 0.));
	visAttMap[2] = new G4VisAttributes(G4Color(1., 0., 0.));
	visAttMap[3] = new G4VisAttributes(G4Color(1., 0., 0.));
	visAttMap[4] = new G4VisAttributes(G4Color(1., 0., 0.));
	visAttMap[5] = new G4VisAttributes(G4Color(0.5, 0.5, 0.5, 0.5));
}

TETParameterisation_EyePlaque::~TETParameterisation_EyePlaque()
{}

G4VSolid* TETParameterisation_EyePlaque::ComputeSolid(
    		       const G4int copyNo, G4VPhysicalVolume* )
{
	// return G4Tet*
	return tetVec[copyNo];
}

void TETParameterisation_EyePlaque::ComputeTransformation(
                   const G4int,G4VPhysicalVolume*) const
{}

G4Material* TETParameterisation_EyePlaque::ComputeMaterial(const G4int copyNo,
                                                 G4VPhysicalVolume* phy,
                                                 const G4VTouchable* )
{
   // set the colour for each organ if visualization is required

	phy->GetLogicalVolume()->SetVisAttributes(visAttMap[idVec[copyNo]]);
	phy->GetLogicalVolume()->SetMaterial(matMap[idVec[copyNo]]);

	// return the material data for each material index
	return matMap[idVec[copyNo]];
}


