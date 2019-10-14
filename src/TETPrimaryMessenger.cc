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
// TETPrimaryMessenger.cc
// \author Haegin Han
//

#include "TETPrimaryGeneratorAction.hh"
#include "TETRunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4RunManager.hh"

TETPrimaryMessenger::TETPrimaryMessenger(TETPrimaryGeneratorAction* _primary)
:G4UImessenger(), fPrimary(_primary), fPrimaryDir(0), fBeamDirCmd(0)
{
	fPrimaryDir = new G4UIdirectory("/beam/");
	fBeamDirCmd = new G4UIcmdWithAString("/beam/dir", this);
	fBeamDirCmd->SetCandidates("AP PA LLAT RLAT ROT ISO");
}

TETPrimaryMessenger::~TETPrimaryMessenger() {
	delete fPrimaryDir;
	delete fBeamDirCmd;
}

void TETPrimaryMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	fPrimary->SetExternalBeam();
	if(command == fBeamDirCmd){
		if(newValue=="AP")	      	fPrimary->SetBeamDirection(AP);
		else if(newValue=="PA")	    fPrimary->SetBeamDirection(PA);
		else if(newValue=="RLAT")	fPrimary->SetBeamDirection(RLAT);
		else if(newValue=="LLAT")	fPrimary->SetBeamDirection(LLAT);
		else if(newValue=="ROT")	fPrimary->SetBeamDirection(ROT);
		else if(newValue=="ISO")	fPrimary->SetBeamDirection(ISO);
	}
}

