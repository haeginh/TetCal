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
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4RunManager.hh"
#include <sstream>
#include <vector>
#include "../include/VoxelRunAction.hh"

TETPrimaryMessenger::TETPrimaryMessenger(TETPrimaryGeneratorAction* _primary)
:G4UImessenger(), fPrimary(_primary)
{
	fExternalDir = new G4UIdirectory("/external/");
	fBeamDirCmd = new G4UIcmdWithAString("/external/dir", this);
	fBeamDirCmd->SetCandidates("AP PA LLAT RLAT ROT ISO");

}

TETPrimaryMessenger::~TETPrimaryMessenger() {
	delete fExternalDir;
	delete fBeamDirCmd;
}

void TETPrimaryMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if(command == fBeamDirCmd){
		fPrimary->SetExternalBeam();
		fPrimary->SetSourceName(newValue);
		ExternalBeam* fExternal = fPrimary->GetExternalBeamGenerator();
		if(newValue=="AP")	      	fExternal->SetBeamDirection(AP);
		else if(newValue=="PA")	    fExternal->SetBeamDirection(PA);
		else if(newValue=="RLAT")	fExternal->SetBeamDirection(RLAT);
		else if(newValue=="LLAT")	fExternal->SetBeamDirection(LLAT);
		else if(newValue=="ROT")	fExternal->SetBeamDirection(ROT);
		else if(newValue=="ISO")	fExternal->SetBeamDirection(ISO);
	}
}

