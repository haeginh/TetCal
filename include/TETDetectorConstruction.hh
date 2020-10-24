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
// TETDetectorConstruction.hh
// \file   MRCP_GEANT4/External/include/TETDetectorConstruction.hh
// \author Haegin Han
//

#ifndef TETDetectorConstruction_h
#define TETDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

#include <cmath>

#include "globals.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tet.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

#include "PSEnergyDeposit.hh"
#include "TETModelImport.hh"
#include "TETParameterisation.hh"

// *********************************************************************
// This is UserDetectorConstruction class that defines geometry
// -- Construct: construct Geometry by three methods listed below.
//  └-- SetupWorldGeometry: Defines the world box (10*10*10 m3) and,
//                          phantom container which has 10 cm-margins from
//                          the bounding box of phantom
//  └-- ConstructPhantom: Define the phantom geometry by using
//                        G4PVParameterised class
//  └-- PrintPhantomInformation: Print overall phantom information
//
// -- ConstructSDandField: Setup the MultiFunctionalDetector with energy
//                         deposition scorer, and attach it to phantom
//                         geometry
// *********************************************************************

class TETDetectorConstruction : public G4VUserDetectorConstruction
{
public:
	TETDetectorConstruction(TETModelImport* tetData);
	virtual ~TETDetectorConstruction();

	virtual G4VPhysicalVolume* Construct();
	virtual void ConstructSDandField();

private:
	void SetupWorldGeometry();
	void ConstructPhantom();
	void PrintPhantomInformation();

	G4VPhysicalVolume* worldPhysical;
	G4LogicalVolume*   container_logic;
    G4VPhysicalVolume* container_phy;

	TETModelImport*    tetData;
	G4ThreeVector      phantomSize;
	G4ThreeVector      phantomBoxMin, phantomBoxMax;
	G4int              nOfTetrahedrons;

	G4LogicalVolume*   tetLogic;

    //4D cal
public:
    bool DeformToBVHFrame(G4int frameNo){
        if(!tetData->Deform(frameNo)) return false;
        G4ThreeVector max = tetData->GetPhantomBoxMax();
        G4ThreeVector min = tetData->GetPhantomBoxMin();
        G4double dimX = max.getX()>-min.getX()? max.getX():-min.getX();
        G4double dimY = max.getY()>-min.getY()? max.getY():-min.getY();
        G4double dimZ = max.getZ()>-min.getZ()? max.getZ():-min.getZ();

        G4Box* phantomBox = dynamic_cast<G4Box*>(container_logic->GetSolid());
        phantomBox->SetXHalfLength(dimX);
        phantomBox->SetYHalfLength(dimY);
        phantomBox->SetZHalfLength(dimZ);
        container_phy->SetTranslation(tetData->GetTranslation(frameNo));
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        currentFrame = frameNo;
    }
    G4int GetCurrentFrameNo(){return currentFrame;}
private:
    G4int currentFrame;
};

#endif
