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

#include "../include/Run.hh"

Run::Run(TETModelImport* tetData)
:G4Run()
{
	fCollID
	= G4SDManager::GetSDMpointer()->GetCollectionID("PhantomSD/eDep");

	organ2dose = tetData->GetDoseMap();

	auto massMap  = tetData->GetMassMap();

	doseOrganized = tetData->DoseWasOrganized();

	//initialize edepMap
    if(!doseOrganized) for(auto itr:massMap) edepMap[itr.first] = {0.,0.};
	else               for(auto itr:organ2dose) edepMap[itr.first] = {0.,0.};
}

Run::~Run()
{
	edepMap.clear();
}

void Run::RecordEvent(const G4Event* event)
{
	// Hits collections
	//
	G4HCofThisEvent* HCE = event->GetHCofThisEvent();
	if(!HCE) return;

	//other doses
	G4THitsMap<G4double>* evtMap =
			static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID));
	auto doseMap = *evtMap->GetMap();
	if(!doseOrganized){
		for(auto itr:doseMap){
			edepMap[itr.first].first += *itr.second;
			edepMap[itr.first].second += (*itr.second)*(*itr.second);
		}
		return;
	}

	//for the organized doses
	std::map<G4int, G4double> edepSum;
	for (auto itr : doseMap) {
		for(auto doseID:organ2dose[itr.first])
			edepSum[doseID]  += *itr.second;
	}
	//organize
	for(auto edep:edepSum){
		edepMap[edep.first].first += edep.second;                 //sum
		edepMap[edep.first].second += edep.second * edep.second;  //square sum
	}
}

void Run::Merge(const G4Run* run)
{
	const Run* localRun = static_cast<const Run*>(run);
	// merge the data from each thread
	EDEPMAP localMap = localRun->edepMap;

	primary = localRun->primary;
	dir = localRun->dir;
	primaryE = localRun->primaryE;
	beamArea = localRun->beamArea;

	for(auto itr : localMap){
		edepMap[itr.first].first  += itr.second.first;
		edepMap[itr.first].second += itr.second.second;
	}

	G4Run::Merge(run);
}

