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
// ActionInitialization.cc
// \file   MRCP_GEANT4/External/src/ActionInitialization.cc
// \author Haegin Han
//

#include "ActionInitialization.hh"
#include "GpsPrimaryGeneratorAction.hh"

ActionInitialization::ActionInitialization(TETModelImport* _tetData, G4String _output, G4Timer* _init, G4bool _useGPS)
 : G4VUserActionInitialization(), tetData(_tetData), output(_output), initTimer(_init), useGPS(_useGPS)
{}

ActionInitialization::~ActionInitialization()
{}

void ActionInitialization::BuildForMaster() const
{
	SetUserAction(new RunAction(tetData, output, initTimer, useGPS));
}

void ActionInitialization::Build() const
{
	// initialise UserAction classes
	if(useGPS) SetUserAction(new GpsPrimaryGeneratorAction());
	else SetUserAction(new PrimaryGeneratorAction(tetData));
	SetUserAction(new RunAction(tetData, output, initTimer, useGPS));
	SetUserAction(new TETSteppingAction);
}  

