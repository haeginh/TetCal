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
// TETRunAction.cc
// \file   MRCP_GEANT4/External/src/TETRunAction.cc
// \author Haegin Han
//

#include "TETRunAction.hh"
#include "G4Timer.hh"
#include <iostream>

TETRunAction::TETRunAction(TETModelImport* _tetData, G4String _output, G4Timer* _init)
:tetData(_tetData), fRun(0), numOfEvent(0), runID(0), outputFile(_output), initTimer(_init), runTimer(0),
 primaryEnergy(-1.), beamArea(-1.), isExternal(true)
{
	if(!isMaster) return;

	runTimer = new G4Timer;
	std::ofstream ofs(outputFile+".node");
	auto skinNodes = tetData->GetSkinNodes();
	ofs<<skinNodes.size()<<"  3  0  0"<<G4endl;
	for(size_t i=0;i<skinNodes.size();i++){
		ofs<<i<<" "<<skinNodes[i].getX()/cm<<" "<<skinNodes[i].getY()/cm<<" "<<skinNodes[i].getZ()/cm<<G4endl;
	}
	ofs.close();
}

TETRunAction::~TETRunAction()
{}

G4Run* TETRunAction::GenerateRun()
{
	// generate run
	fRun = new TETRun(tetData);
	return fRun;
}


void TETRunAction::BeginOfRunAction(const G4Run* aRun)
{
	// print the progress at the interval of 10%
	numOfEvent=aRun->GetNumberOfEventToBeProcessed();
	G4RunManager::GetRunManager()->SetPrintProgress(int(numOfEvent*0.1));
//		    FILE* file = fopen("/proc/self/status", "r");
//		    G4String result;
//		    char line[128];
//
//		    while (fgets(line, 128, file) != NULL){
//		        if (strncmp(line, "VmRSS:", 6) == 0){
//		            result = G4String(line);
//		            break;
//		        }
//		    }
//		    G4cout<<result;
//		    fclose(file);
	if(isMaster){
		initTimer->Stop();
		runTimer->Start();
	}

	const TETPrimaryGeneratorAction* primary =
			dynamic_cast<const TETPrimaryGeneratorAction*>(G4RunManager::GetRunManager()
			->GetUserPrimaryGeneratorAction());
	if(!primary) return;
	primaryParticle = primary->GetParticleGun()->GetParticleDefinition()->GetParticleName();
	primarySourceName = primary->GetSourceName();
	primaryEnergy = primary->GetParticleGun()->GetParticleEnergy();
	beamArea = primary->GetExternalBeamGenerator()->GetBeamArea();
	isExternal = primary-> GetSourceGenerator()->IsExternal();
	fRun->SetPrimary(primaryParticle, primarySourceName, primaryEnergy, beamArea, isExternal);

}

void TETRunAction::EndOfRunAction(const G4Run* aRun)
{
	// print the result only in the Master
	if(!isMaster) return;
	runTimer->Stop();

	// get the run ID
	runID = aRun->GetRunID();

	//get primary info
	primaryParticle = fRun->GetParticleName();
	primarySourceName  = fRun->GetBeamDirName();
	primaryEnergy   = fRun->GetBeamEnergy();
	beamArea        = fRun->GetBeamArea();
	isExternal      = fRun->GetIsExternal();


	// Print the run result by G4cout and std::ofstream
	//

	// set doses
	SetDoses();
/*
	// print by G4cout
	if(isExternal) PrintResultExternal(G4cout);
	else           PrintResultInternal(G4cout);
*/
	// print by std::ofstream
	G4String outName = outputFile + "_" + std::to_string(aRun->GetRunID())+".ele";
	std::ofstream ofs(outName);
	PrintResult(ofs);
	ofs.close();
	G4String errorName = outputFile + "_" + std::to_string(aRun->GetRunID())+".txt";
	std::ofstream ofs2(errorName);
	PrintErrors(ofs2);
	ofs2.close();
/*	if(isExternal)PrintLineExternal(ofs);
	else          PrintLineInternal(ofs);
	ofs.close();*/

	initTimer->Start();
}

void TETRunAction::SetDoses()
{
	doseValues.clear(); doseErrors.clear();
	EDEPMAP edepMap = *fRun->GetEdepMap();

	for(G4int i=0;i<tetData->GetNumSkinTet();i++)
	{
		G4double mass = tetData->GetTetrahedron(tetData->Convert2wholeE(i))->GetCubicVolume()
				        *tetData->GetMaterial(126)->GetDensity();
		G4double meanDose    = edepMap[i].first  / mass / numOfEvent;
		G4double squareDoese = edepMap[i].second / (mass*mass);
		G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
		G4double relativeE   = sqrt(variance)/meanDose;
		doseValues.push_back(meanDose);
		doseErrors.push_back(relativeE);
	}
}
void TETRunAction::PrintResult(std::ostream &out)
{
	out<<tetData->GetNumSkinTet()<<"  4  1"<<G4endl;
	auto eleVec = tetData->GetSkinEle();
	for(G4int i=0;i<tetData->GetNumSkinTet();i++){
		out<<i<<" "<<eleVec[i][0]<<" "<<eleVec[i][1]<<" "<<eleVec[i][2]<<" "<<eleVec[i][3]<<" "<<doseValues[i]/(joule/kg)<<G4endl;
	}
}

void TETRunAction::PrintErrors(std::ostream &out){
	for(G4int i=0;i<tetData->GetNumSkinTet();i++){
		out<<i<<"\t"<<doseValues[i]/(joule/kg)<<"\t"<<doseErrors[i]<<G4endl;
	}
}
