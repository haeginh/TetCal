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
:tetData(_tetData), fRun(0), numOfEvent(0), runID(0), outputFile(_output), initTimer(_init), runTimer(0)
{
	if(isMaster){
		runTimer = new G4Timer;
		std::ofstream ofs(outputFile);
		auto massMap = tetData->GetMassMap();
		ofs<<"[External: pGycm2 / Internal: SAF]"<<G4endl;
		ofs<<"run#\tnps\tinitT\trunT\tparticle\tsource\tenergy[MeV]\t";
		for(auto itr : massMap)
			ofs<<itr.first<<"\t"<<itr.second/g<<"\t";
		ofs<<G4endl;
		ofs.close();
	}
}

TETRunAction::~TETRunAction()
{}

G4Run* TETRunAction::GenerateRun()
{
	// generate run
	fRun = new TETRun();
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
	fRun->SetPrimary(primaryParticle, primarySourceName, primaryEnergy, beamArea);
	isExternal = primary-> GetSourceGenerator() ->IsExternal();
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

	// Print the run result by G4cout and std::ofstream
	//
	// print by G4cout
	PrintResult(G4cout);

	// print by std::ofstream
	std::ofstream ofs(outputFile.c_str(), std::ios::app);
	if(isExternal)PrintLineExternal(ofs);
	else          PrintLineInternal(ofs);
	ofs.close();
	initTimer->Start();
}

void TETRunAction::PrintResult(std::ostream &out)
{
	// Print run result
	//
	using namespace std;
	EDEPMAP edepMap = *fRun->GetEdepMap();

	out << G4endl
	    << "=====================================================================" << G4endl
	    << " Run #" << runID << " / Number of event processed : "<< numOfEvent     << G4endl
	    << "=====================================================================" << G4endl
		<< " Init time: " << initTimer->GetRealElapsed() << " s / Run time: "<< runTimer->GetRealElapsed()<<" s"<< G4endl
	    << "=====================================================================" << G4endl
	    << "organ ID| "
		<< setw(19) << "Organ Mass (g)"
		<< setw(19) << "Dose (Gy*cm2)"
		<< setw(19) << "Relative Error" << G4endl;

	out.precision(3);
	auto massMap = tetData->GetMassMap();
	for(auto itr : massMap){
		G4double meanDose    = edepMap[itr.first].first  / itr.second / numOfEvent;
		G4double squareDoese = edepMap[itr.first].second / (itr.second*itr.second);
		G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
		G4double relativeE   = sqrt(variance)/meanDose;

		out << setw(8)  << itr.first << "| "
			<< setw(19) << fixed      << itr.second/g;
		out	<< setw(19) << scientific << meanDose/(joule/kg) * beamArea/cm2;
		out	<< setw(19) << fixed      << relativeE << G4endl;
	}
	out << "=====================================================================" << G4endl << G4endl;
}

void TETRunAction::PrintLineExternal(std::ostream &out)
{
	// Print run result
	//
	using namespace std;
	EDEPMAP edepMap = *fRun->GetEdepMap();

	out << runID << "\t" <<numOfEvent<<"\t"<< initTimer->GetRealElapsed() << "\t"<< runTimer->GetRealElapsed()<<"\t"
		<< primaryParticle << "\t" <<primarySourceName<< "\t" << primaryEnergy/MeV << "\t";
	auto massMap = tetData->GetMassMap();
	for(auto itr : massMap){
		G4double meanDose    = edepMap[itr.first].first  / itr.second / numOfEvent;
		G4double squareDoese = edepMap[itr.first].second / (itr.second*itr.second);
		G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
		G4double relativeE   = sqrt(variance)/meanDose;

		out << meanDose*1e12/(joule/kg) * beamArea/cm2 <<"\t" << relativeE << "\t";
	}
	out<<G4endl;
}

void TETRunAction::PrintLineInternal(std::ostream &out)
{
	// Print run result
	//
	using namespace std;
	EDEPMAP edepMap = *fRun->GetEdepMap();

	out << runID << "\t" <<numOfEvent<<"\t"<< initTimer->GetRealElapsed() << "\t"<< runTimer->GetRealElapsed()<<"\t"
		<< primaryParticle << "\t" <<primarySourceName<< "\t" << primaryEnergy/MeV << "\t";
	auto massMap = tetData->GetMassMap();
	for(auto itr : massMap){
		G4double meanDose    = edepMap[itr.first].first   / numOfEvent;
		G4double squareDoese = edepMap[itr.first].second;
		G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
		G4double relativeE   = sqrt(variance)/meanDose;

		out << meanDose/primaryEnergy/itr.second /(1./kg) <<"\t" << relativeE << "\t";
	}
	out<<G4endl;
}
