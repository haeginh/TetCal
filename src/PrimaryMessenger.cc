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
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4RunManager.hh"
#include <sstream>
#include <vector>
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"

PrimaryMessenger::PrimaryMessenger(PrimaryGeneratorAction* _primary)
:G4UImessenger(), fPrimary(_primary)
{
	fBeamDir = new G4UIdirectory("/beam/");
	fSpecDirCmd = new G4UIcmdWithAString("/beam/specDir", this);
	fPeakEnergyCmd = new G4UIcmdWithAnInteger("/beam/kVp", this);
	fFanAngleCmd = new G4UIcmdWithADoubleAndUnit("/beam/fanAngle", this);
	fFanAngleCmd->SetDefaultUnit("degree");
	fConeAngleCmd = new G4UIcmdWithADoubleAndUnit("/beam/coneAngle", this);
	fConeAngleCmd->SetDefaultUnit("degree");
	fRadiusCmd = new G4UIcmdWithADoubleAndUnit("/beam/radius", this);
	fRadiusCmd->SetDefaultUnit("cm");
	fLowerBoundCmd = new G4UIcmdWithADoubleAndUnit("/beam/lower", this);
	fLowerBoundCmd->SetDefaultUnit("cm");
	fUpperBoundCmd = new G4UIcmdWithADoubleAndUnit("/beam/upper", this);
	fUpperBoundCmd->SetDefaultUnit("cm");
}

PrimaryMessenger::~PrimaryMessenger() {
	delete fBeamDir;
	delete fSpecDirCmd;
	delete fFanAngleCmd;
	delete fConeAngleCmd;
	delete fRadiusCmd;
	delete fLowerBoundCmd;
	delete fUpperBoundCmd;
}

void PrimaryMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if(command == fSpecDirCmd){
		fPrimary->SetSpecDir(newValue);
	}
	else if(command == fPeakEnergyCmd){
		fPrimary->SetPeakEnergy(fPeakEnergyCmd->GetNewIntValue(newValue));
	}
	else if(command == fFanAngleCmd){
		fPrimary->SetFanAngle(fFanAngleCmd->GetNewDoubleValue(newValue));
	}
	else if(command == fConeAngleCmd){
		fPrimary->SetConeAngle(fConeAngleCmd->GetNewDoubleValue(newValue));
	}
    else if(command == fRadiusCmd){
        fPrimary->SetRadius(fRadiusCmd->GetNewDoubleValue(newValue));
    }
    else if(command == fLowerBoundCmd){
        fPrimary->SetLowerBound(fLowerBoundCmd->GetNewDoubleValue(newValue));
    }
    else if(command == fUpperBoundCmd){
        fPrimary->SetUpperBound(fUpperBoundCmd->GetNewDoubleValue(newValue));
    }
}

