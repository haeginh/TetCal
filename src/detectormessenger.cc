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
/// \file RE06/src/RE06DetectorMessenger.cc
/// \brief Implementation of the RE06DetectorMessenger class
//
//

#include "detectormessenger.hh"

#include "TETDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"


DetectorMessenger::DetectorMessenger(TETDetectorConstruction* det)
 : G4UImessenger(),
   fDetector(det),
   fDirectory(0),
   fDeformCmd(0)
{
    fDirectory = new G4UIdirectory("/4D/");
    fDeformCmd = new G4UIcmdWithAnInteger("/4D/deform", this);
    fDeformCmd->AvailableForStates(G4State_Idle);
}

DetectorMessenger::~DetectorMessenger(){
    delete fDirectory;
    delete fDeformCmd;
}

void DetectorMessenger::SetNewValue(G4UIcommand * command, G4String newValue)
{
    if(command==fDeformCmd) fDetector->DeformToBVHFrame(fDeformCmd->GetNewIntValue(newValue));
}

G4String DetectorMessenger::GetCurrentValue(G4UIcommand *command)
{
    G4String ans;
    if(command==fDeformCmd) ans=fDeformCmd->ConvertToString(fDetector->GetCurrentFrameNo());
    return ans;
}
