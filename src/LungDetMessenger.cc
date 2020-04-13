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
#include "LungDetMessenger.hh"
#include "LungParallelDetCon.hh"

LungDetMessenger::LungDetMessenger(LungParallelDetCon* _lungDet)
:G4UImessenger(), lungDet(_lungDet)
{
    fLungDetDir = new G4UIdirectory("/lung/");
    fVolChkCmd  = new G4UIcmdWithAString("/lung/volchk", this);
    fSamplingCmd= new G4UIcmdWithAnInteger("/lung/sampling", this);
    fBBbasCmd   = new G4UIcmdWithADoubleAndUnit("/lung/BB-bas", this);
    fBBsecCmd   = new G4UIcmdWithADoubleAndUnit("/lung/BB-sec", this);
    fbbsecCmd   = new G4UIcmdWithADoubleAndUnit("/lung/bb-sec", this);
}

LungDetMessenger::~LungDetMessenger() {
    delete fLungDetDir;
    delete fVolChkCmd;
    delete fSamplingCmd;
    delete fBBbasCmd;
    delete fBBsecCmd;
    delete fbbsecCmd;
}

void LungDetMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if(command == fVolChkCmd)
        lungDet->SetVolChkName(newValue);
    else if(command == fSamplingCmd)
        lungDet->SetSamplingNum(fSamplingCmd->GetNewIntValue(newValue));
    else if(command == fBBbasCmd)
        lungDet->SetBBbasVol(fBBbasCmd->GetNewDoubleValue(newValue));
    else if(command == fBBsecCmd)
        lungDet->SetBBsecVol(fBBsecCmd->GetNewDoubleValue(newValue));
    else if(command == fbbsecCmd)
        lungDet->SetbbsecVol(fbbsecCmd->GetNewDoubleValue(newValue));
}

