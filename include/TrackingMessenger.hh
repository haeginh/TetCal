//
// TETTrackingMessenger.hh
//

#ifndef TRACKINGMESSENGER_HH_
#define TRACKINGMESSENGER_HH_ 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAnInteger;
class TETSteppingAction;

class TrackingMessenger: public G4UImessenger
{
public:
	TrackingMessenger(TETSteppingAction* steppingAction);
	virtual ~TrackingMessenger();

	virtual void SetNewValue(G4UIcommand*, G4String);

private:
	TETSteppingAction*         steppingAction;
	G4UIdirectory*             fTrackingDir;
	G4UIcmdWithAnInteger*      fSplitCmd;
};

#endif
