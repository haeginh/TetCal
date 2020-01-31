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

#include "../include/VoxelRunAction.hh"

#include "G4Timer.hh"
#include <iostream>

VoxelRunAction::VoxelRunAction(VOXModelImport* _voxData, G4String _output, G4Timer* _init)
:voxData(_voxData), fRun(0), numOfEvent(0), runID(0), outputFile(_output), initTimer(_init), runTimer(0),
 primaryEnergy(-1.), beamArea(-1.), isExternal(true)
{
	if(!isMaster) return;

	runTimer = new G4Timer;
	std::ofstream ofs(outputFile);

	massMap=voxData->GetMassMap();


	ofs<<"[External: pGycm2 / Internal: SAF (kg-1)]"<<G4endl;
	ofs<<"run#\tnps\tinitT\trunT\tparticle\tsource\tenergy[MeV]\t";
	for(auto itr : massMap)
		ofs<<std::to_string(itr.first)+"_"+voxData->GetOrganName(itr.first)<<"\t"<<itr.second/g<<"\t";
	ofs<<G4endl;
	ofs.close();
}

VoxelRunAction::~VoxelRunAction()
{}

G4Run* VoxelRunAction::GenerateRun()
{
	// generate run
	fRun = new TETRun(voxData);
	return fRun;
}


void VoxelRunAction::BeginOfRunAction(const G4Run* aRun)
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

void VoxelRunAction::EndOfRunAction(const G4Run* aRun)
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

	// print by G4cout
	if(isExternal) PrintResultExternal(G4cout);
	else           PrintResultInternal(G4cout);

	// print by std::ofstream
	std::ofstream ofs(outputFile.c_str(), std::ios::app);
	if(isExternal)PrintLineExternal(ofs);
	else          PrintLineInternal(ofs);
	ofs.close();

	initTimer->Start();
}

void VoxelRunAction::SetDoses()
{
	doseValues.clear(); doseErrors.clear();
	EDEPMAP edepMap = *fRun->GetEdepMap();

	for(auto itr : massMap){
		G4double meanDose    = edepMap[itr.first].first  / itr.second / numOfEvent;
		G4double squareDoese = edepMap[itr.first].second / (itr.second*itr.second);
		G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
		G4double relativeE   = sqrt(variance)/meanDose;
		doseValues.push_back(meanDose);
		doseErrors.push_back(relativeE);
	}
}

void VoxelRunAction::PrintResultExternal(std::ostream &out)
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
	    << setw(27) << "organ ID| "
		<< setw(15) << "Organ Mass (g)"
		<< setw(15) << "Dose (Gy*cm2)"
		<< setw(15) << "Relative Error" << G4endl;

	out.precision(3);

	G4int i=0;

	for(auto itr : massMap){
		out << setw(25) << to_string(itr.first)+"_"+voxData->GetOrganName(itr.first)<< "| ";
		out	<< setw(15) << fixed      << itr.second/g;
		out	<< setw(15) << scientific << doseValues[i]/(joule/kg)*beamArea/cm2;
		out	<< setw(15) << fixed      << doseErrors[i++] << G4endl;
	}

	out << "=====================================================================" << G4endl << G4endl;
}

void VoxelRunAction::PrintResultInternal(std::ostream &out)
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
	    << setw(27) << "organ ID| "
		<< setw(15) << "Organ Mass (g)"
		<< setw(15) << "SAF (kg-1)"
		<< setw(15) << "Relative Error" << G4endl;

	out.precision(3);

	G4int i=0;
	for(auto itr : massMap){
		out << setw(25) <<to_string(itr.first)+"_"+voxData->GetOrganName(itr.first)<< "| ";
		out	<< setw(15) << fixed      << itr.second/g;
		out	<< setw(15) << scientific << doseValues[i]/primaryEnergy/(1./kg);
		out	<< setw(15) << fixed      << doseErrors[i++] << G4endl;
	}

	out << "=====================================================================" << G4endl << G4endl;
}

void VoxelRunAction::PrintLineExternal(std::ostream &out)
{
	// Print run result
	//
	using namespace std;
	EDEPMAP edepMap = *fRun->GetEdepMap();

	out << runID << "\t" <<numOfEvent<<"\t"<< initTimer->GetRealElapsed() << "\t"<< runTimer->GetRealElapsed()<<"\t"
		<< primaryParticle << "\t" <<primarySourceName<< "\t" << primaryEnergy/MeV << "\t";

	for(size_t i=0;i<doseValues.size();i++){
		out << doseValues[i]*1e12/(joule/kg) * beamArea/cm2 <<"\t" << doseErrors[i] << "\t";
	}
	out<<G4endl;
}

void VoxelRunAction::PrintLineInternal(std::ostream &out)
{
	// Print run result
	//
	using namespace std;
	EDEPMAP edepMap = *fRun->GetEdepMap();

	out << runID << "\t" <<numOfEvent<<"\t"<< initTimer->GetRealElapsed() << "\t"<< runTimer->GetRealElapsed()<<"\t"
		<< primaryParticle << "\t" <<primarySourceName<< "\t" << primaryEnergy/MeV << "\t";

	for(size_t i=0;i<doseValues.size();i++){
		out << doseValues[i]/primaryEnergy/(1./kg) <<"\t" << doseErrors[i] << "\t";
	}
	out<<G4endl;
}

