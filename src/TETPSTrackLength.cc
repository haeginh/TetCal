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
// TETPSEnergyDeposit.cc
// \file   MRCP_GEANT4/External/src/TETPSEnergyDeposit.cc
// \author Haegin Han
//

#include "TETPSTrackLength.hh"

TETPSTrackLength::TETPSTrackLength(G4String name, VOXModelImport* _voxData)
  :G4PSTrackLength(name), voxData(_voxData)
{}

TETPSTrackLength::~TETPSTrackLength()
{}

G4int TETPSTrackLength::GetIndex(G4Step* aStep)
{
    G4int iz = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(0);
    G4int ix = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(1);
    G4int iy = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(2);
    G4int idx =  voxData->GetVoxelData(ix,iy,iz);
    if(idx==0) return 0;
    else       return 1;
}
