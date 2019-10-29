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
// TETRun.cc
// \file   MRCP_GEANT4/External/src/TETRun.cc
// \author Haegin Han
//

#include "TETRun.hh"

TETRun::TETRun()
:G4Run()
{
	fCollID
	= G4SDManager::GetSDMpointer()->GetCollectionID("PhantomSD/tLength");
}

TETRun::~TETRun()
{
	lengthBin.clear();
}

void TETRun::RecordEvent(const G4Event* event)
{
	// Hits collections
	//
	G4HCofThisEvent* HCE = event->GetHCofThisEvent();
	if(!HCE) return;

	G4THitsMap<G4double>* evtMap =
			static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID));

	auto doseMap = *evtMap->GetMap();
	lengthBin[floor(*doseMap[0])]++;
}

void TETRun::Merge(const G4Run* run)
{
	const TETRun* localRun = static_cast<const TETRun*>(run);
	// merge the data from each thread
	LENGHBIN localMap = localRun->lengthBin;

	source = localRun->source;
	dir = localRun->dir;
	for(auto itr : localMap){
		lengthBin[itr.first]  += itr.second;
	}

	G4Run::Merge(run);
}






