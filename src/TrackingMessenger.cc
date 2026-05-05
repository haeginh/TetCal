//
// TETTrackingMessenger.cc
//

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "TETSteppingAction.hh"
#include "TrackingMessenger.hh"

TrackingMessenger::TrackingMessenger(TETSteppingAction* _steppingAction)
:G4UImessenger(), steppingAction(_steppingAction), fSplitCmd(0)
{
	fTrackingDir = new G4UIdirectory("/tracking/");
	fSplitCmd = new G4UIcmdWithAnInteger("/tracking/split", this);
	fSplitCmd->SetParameterName("splitValue", false);
	fSplitCmd->AvailableForStates(G4State_Idle, G4State_Init, G4State_PreInit, G4State_GeomClosed);
}

TrackingMessenger::~TrackingMessenger() {
	delete fSplitCmd;
	delete fTrackingDir;
}

void TrackingMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if(command == fSplitCmd){
		G4int splitNum = fSplitCmd->GetNewIntValue(newValue);
		steppingAction->SetSplitNum(splitNum);
		G4cout << "The number of split particles is set to " << splitNum << G4endl;
	}
}
