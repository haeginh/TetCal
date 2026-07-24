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

#include "../include/PSEnergyDeposit.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

PSEnergyDeposit::PSEnergyDeposit(G4String name, TETModelImport* _tetData, G4bool _scoreByTet)
  :G4VPrimitiveScorer(name), tetData(_tetData), scoreByTet(_scoreByTet), HCID(-1), EvtMap(nullptr)
{}

PSEnergyDeposit::~PSEnergyDeposit()
{}

G4bool PSEnergyDeposit::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	G4double edep = aStep->GetTotalEnergyDeposit();
	if(edep == 0.) return FALSE;

	G4double weight = aStep->GetTrack()->GetWeight();
	G4double weightedEdep = edep * weight;
	EvtMap->add(GetIndex(aStep), weightedEdep);
	return TRUE;
}

G4int PSEnergyDeposit::GetIndex(G4Step* aStep)
{
	G4int copyNo = aStep->GetPreStepPoint()->GetTouchable()->GetCopyNumber();
	if(scoreByTet) return copyNo;
	return tetData->GetMaterialIndex(copyNo);
}

void PSEnergyDeposit::Initialize(G4HCofThisEvent* HCE)
{
	EvtMap = new G4THitsMap<G4double>(detector->GetName(), GetName());
	if(HCID < 0) HCID = GetCollectionID(0);
	HCE->AddHitsCollection(HCID,EvtMap);
}

void PSEnergyDeposit::clear()
{
	EvtMap->clear();
}

void PSEnergyDeposit::PrintAll()
{}
