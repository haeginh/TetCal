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
:G4Run(), fCollID_EyeTrackEnter(-1), fCollID_EyeTrackEnterW(-1),
 fCollID_ShieldTrackEnter(-1), fCollID_ShieldTrackEnterW(-1),
 fCollID_ShieldTrackExit(-1), fCollID_ShieldTrackExitW(-1),
 primaryE(-1.), beamArea(-1.), isExternal(false), useSpec(false), useIAEASpec(false),
 shieldTrackEnter(0.), shieldTrackEnterWeighted(0.), shieldTrackExit(0.), shieldTrackExitWeighted(0.)
{
	G4SDManager* sdManager = G4SDManager::GetSDMpointer();
	fCollID = sdManager->GetCollectionID("PhantomSD/eDep");
	fCollID_Hands = sdManager->GetCollectionID("PhantomSD/eDepByTet");
	fCollID_DRF = sdManager->GetCollectionID("PhantomSD/DRF");
	if(sdManager->FindSensitiveDetector("EyeImportanceSD", false)){
		fCollID_EyeTrackEnter = sdManager->GetCollectionID("EyeImportanceSD/TrackEnter");
		fCollID_EyeTrackEnterW = sdManager->GetCollectionID("EyeImportanceSD/TrackEnterW");
	}
	if(sdManager->FindSensitiveDetector("ShieldImportanceSD", false)){
		fCollID_ShieldTrackEnter = sdManager->GetCollectionID("ShieldImportanceSD/TrackEnter");
		fCollID_ShieldTrackEnterW = sdManager->GetCollectionID("ShieldImportanceSD/TrackEnterW");
		fCollID_ShieldTrackExit = sdManager->GetCollectionID("ShieldImportanceSD/TrackExit");
		fCollID_ShieldTrackExitW = sdManager->GetCollectionID("ShieldImportanceSD/TrackExitW");
	}

	organ2dose = tetData->GetDoseMap();
	handTet2dose = tetData->GetHandTetDoseMap();

	auto massMap  = tetData->GetMassMap();
	auto handDoseMassMap = tetData->GetHandDoseMassMap();
	auto rbmRatio = tetData->GetRBMratio();
	auto bsRatio  = tetData->GetBSratio();

	for(auto rbm:rbmRatio)
		rbmFactor[rbm.first] = rbm.second / massMap[rbm.first];
	for(auto bs:bsRatio)
		bsFactor[bs.first] = bs.second / massMap[bs.first];

	doseOrganized = tetData->DoseWasOrganized();

	//initialize edepMap
	edepMap[-1]={0.,0.};
	edepMap[-2]={0.,0.};
	if(!doseOrganized) for(auto itr:massMap) edepMap[itr.first] = {0.,0.};
	else               for(auto itr:organ2dose) edepMap[itr.first] = {0.,0.};
	for(auto itr:handDoseMassMap) edepMap[itr.first] = {0.,0.};
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

	if(fCollID_EyeTrackEnter >= 0){
		G4THitsMap<G4double>* eyeMap =
				static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_EyeTrackEnter));
		if(eyeMap){
			for(auto itr : *eyeMap->GetMap()){
				eyeTrackEnterMap[itr.first] += *itr.second;
			}
		}
	}
	if(fCollID_EyeTrackEnterW >= 0){
		G4THitsMap<G4double>* eyeMapW =
				static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_EyeTrackEnterW));
		if(eyeMapW){
			for(auto itr : *eyeMapW->GetMap()){
				eyeTrackEnterWeightedMap[itr.first] += *itr.second;
			}
		}
	}
	if(fCollID_ShieldTrackEnter >= 0){
		G4THitsMap<G4double>* map =
				static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_ShieldTrackEnter));
		if(map) for(auto itr : *map->GetMap()) shieldTrackEnter += *itr.second;
	}
	if(fCollID_ShieldTrackEnterW >= 0){
		G4THitsMap<G4double>* map =
				static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_ShieldTrackEnterW));
		if(map) for(auto itr : *map->GetMap()) shieldTrackEnterWeighted += *itr.second;
	}
	if(fCollID_ShieldTrackExit >= 0){
		G4THitsMap<G4double>* map =
				static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_ShieldTrackExit));
		if(map) for(auto itr : *map->GetMap()) shieldTrackExit += *itr.second;
	}
	if(fCollID_ShieldTrackExitW >= 0){
		G4THitsMap<G4double>* map =
				static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_ShieldTrackExitW));
		if(map) for(auto itr : *map->GetMap()) shieldTrackExitWeighted += *itr.second;
	}

	//RBM doses
	G4THitsMap<G4double>* evtMap_DRF =
			static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_DRF));
    auto doseMap_DRF = *evtMap_DRF->GetMap();
    for(auto itr:doseMap_DRF){
		edepMap[-4+itr.first].first  += *itr.second;
		edepMap[-4+itr.first].second += (*itr.second)*(*itr.second);
	}

	//other doses
	G4THitsMap<G4double>* evtMap =
			static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID));
	auto doseMap = *evtMap->GetMap();

	if(!handTet2dose.empty()){
		G4THitsMap<G4double>* evtMap_Hands =
				static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_Hands));
		auto handDoseMap = *evtMap_Hands->GetMap();
		std::map<G4int, G4double> handEdepSum;
		for(auto itr : handDoseMap){
			auto doseIDs = handTet2dose.find(itr.first);
			if(doseIDs == handTet2dose.end()) continue;
			for(auto doseID : doseIDs->second) handEdepSum[doseID] += *itr.second;
		}
		for(auto edep : handEdepSum){
			edepMap[edep.first].first += edep.second;
			edepMap[edep.first].second += edep.second * edep.second;
		}
	}

	if(!doseOrganized){
		for(auto itr:doseMap){
			edepMap[itr.first].first += *itr.second;
			edepMap[itr.first].second += (*itr.second)*(*itr.second);
		}
		G4double rbmDose(0.), bsDose(0.);
		for(auto rbm:rbmFactor){
			if(doseMap.find(rbm.first)==doseMap.end()) continue;
			rbmDose += *doseMap[rbm.first] * rbm.second;
		}
		for(auto bs:bsFactor){
			if(doseMap.find(bs.first)==doseMap.end()) continue;
			bsDose += *doseMap[bs.first] * bs.second;
		}
		edepMap[-2].first+=rbmDose; edepMap[-2].second+=rbmDose*rbmDose;
		edepMap[-1].first+=bsDose; edepMap[-1].second+=bsDose*bsDose;
		return;
	}

	//for the organized doses
	std::map<G4int, G4double> edepSum;
	for (auto itr : doseMap) {
		for(auto doseID:organ2dose[itr.first])
			edepSum[doseID]  += *itr.second;
	}
	for(auto rbm:rbmFactor){
		if(doseMap.find(rbm.first)==doseMap.end()) continue;
		edepSum[-2] += *doseMap[rbm.first] * rbm.second;
	}
	for(auto bs:bsFactor){
		if(doseMap.find(bs.first)==doseMap.end()) continue;
		edepSum[-1] += *doseMap[bs.first] * bs.second;
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
	isExternal = localRun->isExternal;
	useSpec = localRun->useSpec;
	useIAEASpec = localRun->useIAEASpec;

	for(auto itr : localMap){
		edepMap[itr.first].first  += itr.second.first;
		edepMap[itr.first].second += itr.second.second;
	}

	for(auto itr : localRun->eyeTrackEnterMap){
		eyeTrackEnterMap[itr.first] += itr.second;
	}
	for(auto itr : localRun->eyeTrackEnterWeightedMap){
		eyeTrackEnterWeightedMap[itr.first] += itr.second;
	}
	shieldTrackEnter += localRun->shieldTrackEnter;
	shieldTrackEnterWeighted += localRun->shieldTrackEnterWeighted;
	shieldTrackExit += localRun->shieldTrackExit;
	shieldTrackExitWeighted += localRun->shieldTrackExitWeighted;

	G4Run::Merge(run);
}
