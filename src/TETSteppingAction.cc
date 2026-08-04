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

#include "ShieldImportanceTuning.hh"
#include "TETModelImport.hh"

#include "G4EventManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4StackManager.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHandle.hh"
#include "G4Track.hh"
#include "G4VUserTrackInformation.hh"
#include "G4VPhysicalVolume.hh"

#include <cmath>

namespace {
class ShieldSplitTrackInformation : public G4VUserTrackInformation
{
  public:
	explicit ShieldSplitTrackInformation(G4bool split = false)
	: shieldSplit(split)
	{}
	~ShieldSplitTrackInformation() override = default;
	void Print() const override {}
	G4bool IsSplit() const { return shieldSplit; }
	void SetSplit() { shieldSplit = true; }

  private:
	G4bool shieldSplit;
};

ShieldSplitTrackInformation* GetShieldTrackInformation(const G4Track* track)
{
	return dynamic_cast<ShieldSplitTrackInformation*>(track->GetUserInformation());
}

G4bool HasShieldSplitInformation(const G4Track* track)
{
	auto information = GetShieldTrackInformation(track);
	return information && information->IsSplit();
}

void MarkShieldSplit(G4Track* track)
{
	auto information = GetShieldTrackInformation(track);
	if(information){
		information->SetSplit();
	}
	else{
		track->SetUserInformation(new ShieldSplitTrackInformation(true));
	}
}
}

TETSteppingAction::TETSteppingAction(TETModelImport* _tetData)
: G4UserSteppingAction(), kCarTolerance(1.0000000000000002e-07), stepCounter(0),
  checkFlag(0), tetData(_tetData)
{}

TETSteppingAction::~TETSteppingAction()
{}

void TETSteppingAction::UserSteppingAction(const G4Step* step)
{
	ApplyShieldSplitting(step);

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
}

G4bool TETSteppingAction::IsShieldMaterial(G4int materialID) const
{
	return materialID == -1 || materialID == -2 || materialID == -7;
}

void TETSteppingAction::ApplyShieldSplitting(const G4Step* step) const
{
	if(!tetData || !ShieldImportanceTuning::IsEnabled()) return;

	G4Track* track = step->GetTrack();
	if(!track || !track->GetParticleDefinition()) return;
	if(track->GetParticleDefinition()->GetParticleName() != ShieldImportanceTuning::GetParticleName()) return;

	const G4StepPoint* pre = step->GetPreStepPoint();
	const G4StepPoint* post = step->GetPostStepPoint();
	if(!pre || !post) return;
	if(post->GetStepStatus() != fGeomBoundary) return;

	G4int preMaterial = 0;
	G4int postMaterial = 0;
	G4bool hasPreTet = false;
	G4bool hasPostTet = false;

	auto preTouchable = pre->GetTouchableHandle();
	if(preTouchable && preTouchable->GetVolume()){
		if(preTouchable->GetVolume()->GetName() == "wholePhantom"){
			G4int copyNo = preTouchable->GetCopyNumber();
			if(copyNo >= 0 && copyNo < tetData->GetNumTetrahedron()){
				preMaterial = tetData->GetMaterialIndex(copyNo);
				hasPreTet = true;
			}
		}
	}

	auto postTouchable = post->GetTouchableHandle();
	if(postTouchable && postTouchable->GetVolume()){
		if(postTouchable->GetVolume()->GetName() == "wholePhantom"){
			G4int copyNo = postTouchable->GetCopyNumber();
			if(copyNo >= 0 && copyNo < tetData->GetNumTetrahedron()){
				postMaterial = tetData->GetMaterialIndex(copyNo);
				hasPostTet = true;
			}
		}
	}

	G4bool preShield = hasPreTet && IsShieldMaterial(preMaterial);
	G4bool postShield = hasPostTet && IsShieldMaterial(postMaterial);

	// The radioactive source generator creates primary tracks (parent ID 0).
	// Restrict splitting to those primaries at their first entry into a shield;
	// interaction-produced secondary gammas are never split here.
	if(!preShield && postShield && track->GetParentID() == 0){
		ShieldImportanceTuning::AddEnter(track->GetWeight());
		if(HasShieldSplitInformation(track)) return;

		G4double splitValue = ShieldImportanceTuning::GetImportance();
		G4int splitCount = std::max(1, (G4int)std::floor(splitValue + 0.5));
		if(splitCount <= 1) return;

		G4StackManager* stack = G4EventManager::GetEventManager()->GetStackManager();
		if(!stack) return;

		const G4double newWeight = track->GetWeight() / splitCount;
		for(G4int i=0; i<splitCount; i++){
			// Replace the parent with complete G4Track copies.  This preserves all
			// kinematic, timing, polarization, vertex, and touchable state equally
			// for every branch instead of treating the original asymmetrically.
			G4Track* clone = new G4Track(*track);
			clone->SetWeight(newWeight);
			clone->SetTrackID(0);
			clone->SetParentID(track->GetParentID());
			clone->SetTrackStatus(fAlive);
			MarkShieldSplit(clone);
			stack->PushOneTrack(clone);
		}
		track->SetTrackStatus(fStopAndKill);
	}
	else if(preShield && !postShield && HasShieldSplitInformation(track)){
		ShieldImportanceTuning::AddExit(track->GetWeight());
	}
}
