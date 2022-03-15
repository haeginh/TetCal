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
	fCollID_DRF
	= G4SDManager::GetSDMpointer()->GetCollectionID("PhantomSD/DRF");

	organ2dose = tetData->GetDoseMap();

	auto massMap  = tetData->GetMassMap();
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

	source = localRun->source;
	factor = localRun->factor;
	
	for(auto itr : localMap){
		edepMap[itr.first].first  += itr.second.first;
		edepMap[itr.first].second += itr.second.second;
	}

	G4Run::Merge(run);
}

