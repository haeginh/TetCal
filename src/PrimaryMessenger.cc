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
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4RunManager.hh"
#include <sstream>
#include <vector>
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"

PrimaryMessenger::PrimaryMessenger(PrimaryGeneratorAction* _primary)
:G4UImessenger(), fPrimary(_primary)
{
    fExternalDir = new G4UIdirectory("/beam/");
    fBeamDirCmd = new G4UIcmdWithAString("/beam/dir", this);
	fBeamDirCmd->SetCandidates("AP PA LLAT RLAT ROT ISO");

    fBeamSizeCmd = new G4UIcmdWith3VectorAndUnit("/beam/size", this);
    fBeamSizeCmd->SetDefaultUnit("cm");
    fBeamSizeCmd->SetUnitCandidates("cm mm m");

    fBeamCenterCmd = new G4UIcmdWith3VectorAndUnit("/beam/center", this);
    fBeamCenterCmd->SetDefaultUnit("cm");
    fBeamCenterCmd->SetUnitCandidates("cm mm m");
    fBeamCenterCmd->SetDefaultValue(G4ThreeVector());

    fBeamRadiusCmd = new G4UIcmdWithADoubleAndUnit("/beam/radius", this);
    fBeamRadiusCmd->SetDefaultUnit("cm");
    fBeamRadiusCmd->SetUnitCandidates("cm mm m");

    fBeamDefaultCmd = new G4UIcmdWithoutParameter("/beam/default", this);
}

PrimaryMessenger::~PrimaryMessenger() {
	delete fExternalDir;
	delete fBeamDirCmd;
    delete fBeamSizeCmd;
    delete fBeamCenterCmd;
    delete fBeamDefaultCmd;
}

void PrimaryMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if(command == fBeamDirCmd){
        fPrimary->SetSourceName(newValue);
        ExternalBeam* fExternal = fPrimary->GetSourceGenerator();
		if(newValue=="AP")	      	fExternal->SetBeamDirection(AP);
		else if(newValue=="PA")	    fExternal->SetBeamDirection(PA);
		else if(newValue=="RLAT")	fExternal->SetBeamDirection(RLAT);
		else if(newValue=="LLAT")	fExternal->SetBeamDirection(LLAT);
		else if(newValue=="ROT")	fExternal->SetBeamDirection(ROT);
		else if(newValue=="ISO")	fExternal->SetBeamDirection(ISO);
	}
    else if(command == fBeamSizeCmd){
        fPrimary->GetSourceGenerator()->SetBeamXYZ(fBeamSizeCmd->GetNew3VectorValue(newValue));
    }
    else if(command == fBeamCenterCmd){
        fPrimary->GetSourceGenerator()->SetBeamCenter(fBeamCenterCmd->GetNew3VectorValue(newValue));
    }
    else if(command == fBeamRadiusCmd){
        fPrimary->GetSourceGenerator()->SetBeamRadius(fBeamRadiusCmd->GetNewDoubleValue(newValue));
    }
    else if(command == fBeamDefaultCmd){
        fPrimary->GetSourceGenerator()->SetDefaultSize();
    }
}

