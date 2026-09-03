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
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"

PrimaryMessenger::PrimaryMessenger(PrimaryGeneratorAction* _primary)
:G4UImessenger(), fPrimary(_primary)
{
	fExternalDir = new G4UIdirectory("/external/");
	fBeamDirCmd = new G4UIcmdWithAString("/external/dir", this);
	fBeamDirCmd->SetCandidates("AP PA LLAT RLAT ROT ISO");

	fInternalDir      = new G4UIdirectory("/internal/");
	fSourceOrganCmd   = new G4UIcmdWithAString("/internal/source", this);
	fSourceOrganFractionCmd = new G4UIcmdWithAString("/internal/fraction", this);
	fSurfaceSourceCmd = new G4UIcmdWithAString("/internal/surface", this);
	fSourceOrganCmd->SetGuidance("Set the internal source organ IDs.");
	fSourceOrganFractionCmd->SetGuidance(
		"Set one non-negative source weight per /internal/source organ ID.");

	fSpectrumDir = new G4UIdirectory("/spec/");
	fSpectrumSourceCmd = new G4UIcmdWithAString("/spec/input", this);
	fRadCodesCmd = new G4UIcmdWithAString("/spec/RADcodes", this);
}

PrimaryMessenger::~PrimaryMessenger() {
	delete fBeamDirCmd;
	delete fSourceOrganCmd;
	delete fSourceOrganFractionCmd;
	delete fSurfaceSourceCmd;
	delete fSpectrumSourceCmd;
	delete fRadCodesCmd;
	delete fExternalDir;
	delete fInternalDir;
	delete fSpectrumDir;
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
	else if(command == fSourceOrganCmd){
		fPrimary->SetInternalBeam();
		InternalSource* fInternal = fPrimary->GetInternalBeamGenerator();
        if(newValue.substr(0, 1)=="\"") newValue = newValue.substr(1, newValue.size()-2);

		std::stringstream ss(newValue);
		std::vector<G4int> organIDs;
		G4int intTemp;
		while(ss>>intTemp) organIDs.push_back(intTemp);
		if(!ss.eof()){
			G4Exception("PrimaryMessenger::SetNewValue", "InvalidSourceList",
					FatalErrorInArgument, "Invalid value in /internal/source");
			return;
		}
		if(organIDs.empty()){
			G4Exception("PrimaryMessenger::SetNewValue", "EmptySourceList",
					FatalErrorInArgument, "/internal/source requires at least one organ ID");
			return;
		}
		fInternal->SetSource(organIDs);
		fCurrentSourceText = newValue;
		if(fPendingFractions.empty()){
			fPrimary->SetSourceName("(V) "+newValue);
		} else {
			fInternal->SetFractions(fPendingFractions);
			fPendingFractions.clear();
			fPrimary->SetSourceName("(VF) "+newValue);
		}
	}
	else if(command == fSourceOrganFractionCmd){
		fPrimary->SetInternalBeam();
		if(newValue.substr(0, 1)=="\"") newValue = newValue.substr(1, newValue.size()-2);

		std::stringstream ss(newValue);
		std::vector<G4double> fractions;
		G4double doubleTemp;
		while(ss>>doubleTemp) fractions.push_back(doubleTemp);
		if(!ss.eof()){
			G4Exception("PrimaryMessenger::SetNewValue", "InvalidFractionList",
					FatalErrorInArgument, "Invalid value in /internal/fraction");
			return;
		}
		if(fractions.empty()){
			G4Exception("PrimaryMessenger::SetNewValue", "EmptyFractionList",
					FatalErrorInArgument, "/internal/fraction requires at least one weight");
			return;
		}

		InternalSource* fInternal = fPrimary->GetInternalBeamGenerator();
		if(fInternal->HasSource()){
			fInternal->SetFractions(fractions);
			fPrimary->SetSourceName("(VF) "+fCurrentSourceText);
		} else {
			fPendingFractions = fractions;
			G4cout << "Stored " << fractions.size()
			       << " source fractions; waiting for /internal/source" << G4endl;
		}
	}
    else if(command == fSurfaceSourceCmd){
        fPrimary->SetSurfaceSource();
        SurfaceSource* fSurface = fPrimary->GetSurfaceSourceGenerator();
        if(newValue.substr(0, 1)=="\"") newValue = newValue.substr(1, newValue.size()-2);

        fPrimary->SetSourceName("(S) "+newValue);

        std::stringstream ss(newValue);
        std::vector<G4int> organIDs;
        G4int intTemp;
        while(ss>>intTemp) organIDs.push_back(intTemp);
        fSurface->SetSource(organIDs);
    }
	else if(command == fSpectrumSourceCmd){
		fPrimary->SetSpectrumSource(newValue);
	}
	else if(command == fRadCodesCmd){
		fPrimary->SetRadCodes(newValue);
	}
}
