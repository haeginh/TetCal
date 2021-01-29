#ifndef TETDATAMESSENGER_HH
#define TETDATAMESSENGER_HH

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class LungParallelDetCon;

class LungDetMessenger: public G4UImessenger
{
public:
    LungDetMessenger(LungParallelDetCon* _lungDet);
    virtual ~LungDetMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);

private:
    LungParallelDetCon*        lungDet;
    G4UIdirectory*             fLungDetDir;
    G4UIcmdWithAString*        fVolChkCmd;
    G4UIcmdWithADoubleAndUnit* fBBbasCmd;
    G4UIcmdWithADoubleAndUnit* fBBsecCmd;
    G4UIcmdWithADoubleAndUnit* fbbsecCmd;
};


#endif // TETDATAMESSENGER_HH
