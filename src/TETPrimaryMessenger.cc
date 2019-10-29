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
#include <sstream>
#include <vector>

TETPrimaryMessenger::TETPrimaryMessenger(TETPrimaryGeneratorAction* _primary)
:G4UImessenger(), fPrimary(_primary)
{
	fInternalDir      = new G4UIdirectory("/odd/");
	fSourceOrganCmd   = new G4UIcmdWithAString("/odd/organ", this);
	fBeamDirCmd       = new G4UIcmdWithAString("/odd/dir", this);
	fBeamDirCmd->SetCandidates("front back left right top bottom");
}

TETPrimaryMessenger::~TETPrimaryMessenger() {
	delete fInternalDir;
	delete fSourceOrganCmd;
}

void TETPrimaryMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
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
	if(command == fBeamDirCmd){
		if(newValue == "front")
			fPrimary->GetParticleGun()->SetParticleMomentumDirection(G4ThreeVector(0., -1., 0.));
		if(newValue == "back")
			fPrimary->GetParticleGun()->SetParticleMomentumDirection(G4ThreeVector(0., 1., 0.));
		if(newValue == "left")
			fPrimary->GetParticleGun()->SetParticleMomentumDirection(G4ThreeVector(1., 0., 0.));
		if(newValue == "right")
			fPrimary->GetParticleGun()->SetParticleMomentumDirection(G4ThreeVector(-1., 0., 0.));
		if(newValue == "top")
			fPrimary->GetParticleGun()->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
		if(newValue == "bottom")
			fPrimary->GetParticleGun()->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
	}
}

