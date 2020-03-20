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

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4RunManager.hh"
#include <sstream>
#include <vector>
#include "../include/PrimaryGeneratorAction.hh"
#include "../include/RunAction.hh"

PrimaryMessenger::PrimaryMessenger(PrimaryGeneratorAction* _primary)
:G4UImessenger(), fPrimary(_primary)
{
	fExternalDir = new G4UIdirectory("/external/");
	fBeamDirCmd = new G4UIcmdWithAString("/external/dir", this);
	fBeamDirCmd->SetCandidates("AP PA LLAT RLAT ROT ISO");

	fInternalDir      = new G4UIdirectory("/internal/");
	fSourceOrganCmd   = new G4UIcmdWithAString("/internal/source", this);
	fSurfaceSourceCmd = new G4UIcmdWithAString("/internal/surface", this);
}

PrimaryMessenger::~PrimaryMessenger() {
	delete fExternalDir;
	delete fBeamDirCmd;
	delete fInternalDir;
	delete fSourceOrganCmd;
}

void PrimaryMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
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
	if(command == fSourceOrganCmd){
		fPrimary->SetInternalBeam();
		InternalSource* fInternal = (InternalSource*) fPrimary->GetInternalBeamGenerator();
        if(newValue.substr(0, 1)=="\"") newValue = newValue.substr(1, newValue.size()-2);

		fPrimary->SetSourceName("(V) "+newValue);

		std::stringstream ss(newValue);
		std::vector<G4int> organIDs;
		G4int intTemp;
		while(ss>>intTemp) organIDs.push_back(intTemp);
		fInternal->SetSource(organIDs);
	}
    if(command == fSurfaceSourceCmd){
        fPrimary->SetSurfaceSource();
        SurfaceSource* fSurface = (SurfaceSource*) fPrimary->GetSurfaceSourceGenerator();
        if(newValue.substr(0, 1)=="\"") newValue = newValue.substr(1, newValue.size()-2);

        fPrimary->SetSourceName("(S) "+newValue);

        std::stringstream ss(newValue);
        std::vector<G4int> organIDs;
        G4int intTemp;
        while(ss>>intTemp) organIDs.push_back(intTemp);
        fSurface->SetSource(organIDs);
    }
}

