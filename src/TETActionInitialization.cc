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
// TETActionInitialization.cc
// \file   MRCP_GEANT4/External/src/TETActionInitialization.cc
// \author Haegin Han
//

#include "TETActionInitialization.hh"

TETActionInitialization::TETActionInitialization(TETModelImport* _tetData, G4String _output, G4Timer* _init)
 : G4VUserActionInitialization(), tetData(_tetData), output(_output), initTimer(_init)
{}

TETActionInitialization::~TETActionInitialization()
{}

void TETActionInitialization::BuildForMaster() const
{
	SetUserAction(new TETRunAction(tetData, output, initTimer));
}

void TETActionInitialization::Build() const
{
	// initialise UserAction classes
	SetUserAction(new TETPrimaryGeneratorAction());
	SetUserAction(new TETRunAction(tetData, output, initTimer));
	SetUserAction(new TETSteppingAction);
}  

