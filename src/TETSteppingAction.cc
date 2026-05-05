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
// TETSteppingAction.cc
// \file   MRCP_GEANT4/External/src/TETSteppingAction.cc
// \author Haegin Han
//

#include "TETSteppingAction.hh"
#include "G4Gamma.hh"
#include "G4EventManager.hh"
#include "G4StackManager.hh"

TETSteppingAction::TETSteppingAction(TETModelImport* _tetmodel, G4int _splitMaterialID)
: G4UserSteppingAction(), tetmodel(_tetmodel), trackingMessenger(0), splitMaterialID(_splitMaterialID),
  splitNum(1), splitNum_inv(1.), kCarTolerance(1.0000000000000002e-07), stepCounter(0), checkFlag(0)
{
	trackingMessenger = new TrackingMessenger(this);
}

TETSteppingAction::~TETSteppingAction()
{
	delete trackingMessenger;
}

void TETSteppingAction::SetSplitNum(G4int num)
{
	if(num <= 0){
		G4Exception("TETSteppingAction::SetSplitNum","",FatalErrorInArgument,
				G4String("      The number of split particles must be positive: "
						+ std::to_string(num)).c_str());
	}
	splitNum = num;
	splitNum_inv = 1. / splitNum;
}

void TETSteppingAction::UserSteppingAction(const G4Step* step)
{
	// Slightly move the particle when the step length of five continuous steps is
	// shorter than the tolerance (0.1 nm)
	//
	G4Track* theTrack = step->GetTrack();
	G4bool CheckingLength = (step->GetStepLength() < kCarTolerance);
	if(CheckingLength)
	{
		++stepCounter;
		if( checkFlag && stepCounter>=5 )
		{
			// kill the track if the particle is stuck even after the slight move
			// (this hardly occurs)
			theTrack->SetTrackStatus(fStopAndKill);
			stepCounter=0;
			checkFlag=0;
		}
		else if( stepCounter>=5 )
		{
			// if a particle is at the same position (step length < 0.1 nm) for five consecutive steps,
			// slightly move (0.1 nm) the stuck particle in the direction of momentum
			theTrack->SetPosition(theTrack->GetPosition() + theTrack->GetMomentumDirection()*kCarTolerance);
			checkFlag=1;
		}
	}
	else stepCounter=0;

	if(splitMaterialID < 0) return;

	// Remove all secondary particles generated in the selected apron material.
	if (step->GetPostStepPoint()->GetTouchable()->GetHistoryDepth()==2&&
	    tetmodel->GetMaterialIndex(step->GetPostStepPoint()->GetTouchable()->GetCopyNumber()) == splitMaterialID) {
		auto secondaries = step->GetSecondaryInCurrentStep();
		for (auto* secondary : *secondaries) {
			const_cast<G4Track*>(secondary)->SetTrackStatus(fStopAndKill);
		}
	}

	if(splitNum==1) return;
	if(theTrack->GetKineticEnergy() < 30.*keV) return;

	auto* info = dynamic_cast<MyTrackInfo*>(theTrack->GetUserInformation());
	if (info) return;
	if (theTrack->GetParticleDefinition()!=G4Gamma::Gamma()) return;
	if(step->GetPostStepPoint()->GetStepStatus()!=fGeomBoundary) return;

	if(step->GetPreStepPoint()->GetTouchable()->GetCopyNumber()==-1
	  &&tetmodel->GetMaterialIndex(step->GetPostStepPoint()->GetTouchable()->GetCopyNumber())==splitMaterialID){
		for (G4int i = 0; i < splitNum; ++i) {
			G4Track* newTrack = new G4Track(*theTrack);
			newTrack->SetUserInformation(new MyTrackInfo());
			newTrack->SetWeight(theTrack->GetWeight() * splitNum_inv);
			newTrack->SetTrackStatus(fAlive);
			G4EventManager::GetEventManager()->GetStackManager()->PushOneTrack(newTrack);
		}
		theTrack->SetTrackStatus(fStopAndKill);
	}
}
