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
/// \file runAndEvent/TET/src/TETDRFDose.cc
/// \brief Implementation of the TETDRFDose class
//
//
// TETDRFDose
#include "TETDRFDose.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "TETDetectorConstruction.hh"
#include "G4RunManager.hh"
//using namespace std;
///////////////////////////////////////////////////////////////////////////////
// (Description) DRF scorer for adult MRCPs
//
// Created: 2015-12-14  Yeon Soo Yeom.
// Edited : 2019-10-31  Haegin Han.
///////////////////////////////////////////////////////////////////////////////

TETDRFDose::TETDRFDose(G4String name,TETModelImport* _PhantomData)
  :G4VPrimitiveScorer(name), phantomData(_PhantomData), HCID(-1), EvtMap(0)
{
	energyBin={0.01 ,0.015 ,0.02 ,0.03 ,0.04 ,0.05 ,0.06 ,0.08
				,0.1 ,0.15 ,0.2 ,0.3 ,0.4 ,0.5 ,0.6 ,0.8
				,1 ,1.5 ,2 ,3 ,4 ,5 ,6 ,8 ,10};
	gamma = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
	rbmDRF = phantomData->GetRBMDRF();
	bsDRF = phantomData->GetBSDRF();
	rbmRatio = phantomData->GetRBMmap();
	bsRatio = phantomData->GetBSmap();
}

TETDRFDose::~TETDRFDose()
{;}

G4int TETDRFDose::GetIndex(G4Step* aStep)
{
	// return the organ ID (= material index)
	G4int copyNo = aStep->GetPreStepPoint()->GetTouchable()->GetCopyNumber();
	return phantomData->GetMaterialIndex(copyNo);
}

G4bool TETDRFDose::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	if (aStep -> GetTrack() -> GetDefinition() != gamma) return FALSE;
	G4double stepLength = aStep->GetStepLength();
	if (stepLength==0.) return FALSE;
	G4int index = GetIndex(aStep);
	if (rbmRatio.find(index)==rbmRatio.end()) return FALSE;
	G4double CellFlux = stepLength / phantomData->GetVolume(index);
	G4double energy=aStep->GetPreStepPoint()->GetKineticEnergy();
	G4double RBMdose = GetRBMdose(energy, CellFlux, index);
	G4double BSdose = GetBSdose(energy, CellFlux, index);

	if(std::isnan(RBMdose)) RBMdose=0;
	if(std::isnan(BSdose)) BSdose=0;

	EvtMap->add(0, RBMdose);
	EvtMap->add(1, BSdose);

	return TRUE;
}

void TETDRFDose::Initialize(G4HCofThisEvent* HCE)
{
	EvtMap = new G4THitsMap<G4double>(detector->GetName(), GetName());
    if(HCID < 0) HCID = GetCollectionID(0);
    HCE->AddHitsCollection(HCID,EvtMap);
}

void TETDRFDose::EndOfEvent(G4HCofThisEvent*)
{;}

void TETDRFDose::clear()
{
	EvtMap->clear();
}

G4double TETDRFDose::GetRBMdose(G4double energy, G4double cellFlux, G4int organID)
{
	G4int eIdx = FindIndexfromEnergyBin(energy);
	G4double RBM = rbmRatio[organID];
	G4double RBMDRF = rbmDRF[organID][eIdx];
	G4double NextRBMDRF = rbmDRF[organID][eIdx+1];

    if (RBMDRF==0) return 0.0;

    G4double DRF = log10( RBMDRF )
					+ (log10(energy/(energyBin)[eIdx]))
					* log10(NextRBMDRF/RBMDRF)
					/ log10((energyBin)[eIdx+1]/(energyBin)[eIdx]);

	G4double RBMDose = exp10(DRF) * 1e6 * cellFlux * RBM; //Convert to Gy

    return RBMDose;
}
G4double TETDRFDose::GetBSdose(G4double energy, G4double cellFlux, G4int organID){

	G4int eIdx = FindIndexfromEnergyBin(energy);
	G4double BS = bsRatio[organID];
	G4double BSDRF = bsDRF[organID][eIdx];
	G4double NextBSDRF = bsDRF[organID][eIdx+1];

    if (BSDRF==0) return 0.0;

    G4double DRF = log10( BSDRF )
					+ (log10(energy/(energyBin)[eIdx]))
					* log10(NextBSDRF/BSDRF)
					/ log10((energyBin)[eIdx+1]/(energyBin)[eIdx]);

	G4double BSDose = exp10(DRF) * 1e6 * cellFlux * BS; //Convert to Gy

    return BSDose;
}

G4int TETDRFDose::FindIndexfromEnergyBin(G4double energy){
	for(G4int i=0;i<(G4int)energyBin.size()-1;i++) {
		if((energyBin)[i+1] >= energy) return i;
	}
	return (G4int)(energyBin).size()-2;
}


