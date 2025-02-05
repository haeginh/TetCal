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

#include "DosimeterScorer.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "TETDetectorConstruction.hh"
#include "G4RunManager.hh"

//Unit: Gy
DosimeterScorer::DosimeterScorer(G4String name,TETModelImport* _PhantomData)
  :G4VPrimitiveScorer(name), PhantomData(_PhantomData), HCID(-1), HCID1(-1), EvtMap(0)
{
	// energyBin={0.01 ,0.015 ,0.02 ,0.03 ,0.04 ,0.05 ,0.06 ,0.08
	// 			,0.1 ,0.15 ,0.2 ,0.3 ,0.4 ,0.5 ,0.6 ,0.8
	// 			,1 ,1.5 ,2 ,3 ,4 ,5 ,6 ,8 ,10};
	cosID[cos(15*deg)] = 1;
	cosID[cos(30*deg)] = 2;
	cosID[cos(45*deg)] = 3;
	cosID[cos(60*deg)] = 4;
	cosID[cos(75*deg)] = 5;
	std::ifstream ifs("Hp10_DC.dat");
	G4String line;
	std::getline(ifs, line);
	std::stringstream ss(line);
	G4double val;
	std::vector<G4double> bin;
	while(ss>>val) bin.push_back(val*MeV);
	for(G4int i=0;i<cosID.size()+1;i++){
		for(G4int b=0;b<bin.size();b++){
			ifs>>val;
			coeff[bin[b]].push_back(val);
		}
	}
	ifs.close();
	for(G4int b=0;b<bin.size();b++){
		eBinID[bin[b]] = b;
	}

	auto dosiEle = PhantomData->GetDosiEle();
	for(G4int i=0;i<dosiEle.size();i++){
		for(G4int e:dosiEle[i]) dosimeter[e].push_back(i);
	}
	// RBMratio = PhantomData->GetRBMratio();
	// BSratio = PhantomData->GetBSratio();
	gamma = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
}

DosimeterScorer::~DosimeterScorer()
{;}

// G4int DosimeterScorer::GetIndex(G4Step* aStep)
// {
// 	G4int copyNo = aStep->GetPreStepPoint()->GetTouchable()->GetCopyNumber();
// 	return PhantomData->GetMaterialIndex(copyNo);
// }

G4bool DosimeterScorer::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	G4int index = aStep -> GetPostStepPoint() ->GetTouchable() -> GetCopyNumber();
	// std::cout<<index<<std::endl;
	// if(RBMratio.find(index)==RBMratio.end()) return FALSE;
	if (dosimeter.find(index)==dosimeter.end()) return FALSE;
	if (aStep -> GetPreStepPoint() ->GetTouchable() -> GetCopyNumber() >= 0) return FALSE;
	if (aStep -> GetTrack() -> GetDefinition()!=gamma) return FALSE;
	
	G4double cos = aStep->GetTrack()->GetMomentumDirection().dot(PhantomData->GetDosiNorm(index));
	G4double factor = GetFactor(aStep->GetTrack()->GetKineticEnergy(), cos);	
	G4double dose = factor*cos;
	// G4cout<<index<<G4endl;
	for(G4int d:dosimeter[index]){
		// G4cout<<"d"<<d<<":"<<factor<<" "<<cos<< " "<<" "<<dose<<G4endl;
		EvtMap->add(d, dose);
		// G4cout<<"dosimeter: "<<d<<G4endl;
		G4int binID = -(eBinID[coeff.lower_bound(aStep->GetTrack()->GetKineticEnergy())->first]+100*d+1);
		EvtMap->add(binID, cos);
	}

//	G4int copyNo = aStep->GetTrack()->GetOriginTouchable()->Get;
//	G4cout<<PhantomData->GetMaterialIndex(copyNo)<<"\t"<<index<<G4endl;
//	if(PhantomData->GetMaterialIndex(copyNo)!=index) {
//		EvtMap->add(-2, RBMdose);
//		EvtMap->add(-1, BSdose);
//	}

	return TRUE;
}

void DosimeterScorer::Initialize(G4HCofThisEvent* HCE)
{
	EvtMap = new G4THitsMap<G4double>(detector->GetName(), GetName());
	if(HCID < 0) HCID = GetCollectionID(0);
    HCE->AddHitsCollection(HCID,EvtMap);
}

void DosimeterScorer::EndOfEvent(G4HCofThisEvent*)
{;}

void DosimeterScorer::clear()
{
	EvtMap->clear();
}

G4double DosimeterScorer::GetFactor(G4double energy, G4double cos){
	// G4cout<<energy<<" "<<cos<<G4endl;
	auto upperCosIter = cosID.upper_bound(cos); // up (selfX)
	auto upperEnergyIter = coeff.upper_bound(energy); // up (selfX)
	G4double upperE = upperEnergyIter->first;
	G4double upperE_val = upperEnergyIter->second[0];
	G4double upperCos(1), lowerCos(0);
	if(upperEnergyIter==coeff.begin()) { //energy is lower than the smallest energy bin (same as the smallest E)
		G4double upperCos_val(1), lowerCos_val(0);
		if(upperCosIter!=cosID.begin()){lowerCos = prev(upperCosIter)->first; upperCos_val = 1;}
		else if(upperCosIter!=cosID.end()){upperCos = upperCosIter->first; lowerCos_val = 0;}
		G4double angleF = (lowerCos_val+(upperCos_val-lowerCos_val)*(cos-lowerCos)/(upperCos-lowerCos));
		G4double eCoeff = upperE_val*energy/upperE;
		return angleF*eCoeff;
	}
	G4double lowerE = prev(upperEnergyIter)->first;
	G4double lowerE_val = prev(upperEnergyIter)->second[0];
	G4double eCoeff = (lowerE_val+(upperE_val-lowerE_val)*(energy-lowerE)/(upperE-lowerE));
	G4double val_uu, val_ll, val_ul, val_lu; //angle factor (energy, angle)
	if(upperCosIter==cosID.begin()){
		upperCos = upperCosIter->first;
		val_uu = upperEnergyIter->second[upperCosIter->second];
		val_lu = prev(upperEnergyIter)->second[upperCosIter->second];
		G4double angleF_0u = val_lu+(val_uu-val_lu)*(energy-lowerE)/(upperE-lowerE);
		G4double angleF = angleF_0u*cos/upperCos;
		// G4cout<<"1 "<<val_uu<<" "<<val_lu<<" "<<angleF<<" "<<eCoeff<<G4endl;
		return angleF*eCoeff;
	}
	if(upperCosIter==cosID.end()){
		lowerCos = prev(upperCosIter)->first;
		val_ul = upperEnergyIter->second[prev(upperCosIter)->second];
		val_ll = prev(upperEnergyIter)->second[prev(upperCosIter)->second];
		G4double angleF_0l = val_ll+(val_ul-val_ll)*(energy-lowerE)/(upperE-lowerE);
		G4double angleF = angleF_0l+(1-angleF_0l)*(cos-lowerCos)/(1-lowerCos);
		// G4cout<<"2 "<<val_ul<<" "<<val_ll<<" "<<angleF<<" "<<eCoeff<<G4endl;
		return angleF*eCoeff;
	}
	upperCos = upperCosIter->first;
	lowerCos = prev(upperCosIter)->first;
	val_uu = upperEnergyIter->second[upperCosIter->second];
	val_lu = prev(upperEnergyIter)->second[upperCosIter->second];
	G4double angleF_0u = val_lu+(val_uu-val_lu)*(energy-lowerE)/(upperE-lowerE);
	val_ul = upperEnergyIter->second[prev(upperCosIter)->second];
	val_ll = prev(upperEnergyIter)->second[prev(upperCosIter)->second];
	G4double angleF_0l = val_ll+(val_ul-val_ll)*(energy-lowerE)/(upperE-lowerE);
	G4double angleF = angleF_0l+(angleF_0u-angleF_0l)*(cos-lowerCos)/(upperCos-lowerCos);
	// G4cout<<"3 "<<val_uu<<" "<<val_lu<<" "<<val_ul<<" "<<val_ll<<" "<<angleF<<" "<<eCoeff<<G4endl;
	return angleF*eCoeff;
}
