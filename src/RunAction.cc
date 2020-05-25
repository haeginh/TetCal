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

#include "G4Timer.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4Neutron.hh"
#include <iostream>
#include "RunAction.hh"

RunAction::RunAction(TETModelImport* _tetData, G4String _output, G4Timer* _init)
:tetData(_tetData), fRun(0), numOfEvent(0), runID(0), outputFile(_output), initTimer(_init), runTimer(0),
 primaryEnergy(-1.), beamArea(-1.)
{
	if(!isMaster) return;

	runTimer = new G4Timer;
	std::ofstream ofs(outputFile);

	if(tetData->DoseWasOrganized()){
		massMap = tetData->GetDoseMassMap();
		for(auto itr:massMap) nameMap[itr.first] = tetData->GetDoseName(itr.first);
	}
	else{
		massMap = tetData->GetMassMap();
		for(auto itr:massMap) nameMap[itr.first] = tetData->GetMaterial(itr.first)->GetName();
	}

	ofs<<"[External: pGycm2 / Internal: SAF (kg-1)]"<<G4endl;
	ofs<<"run#\tnps\tinitT\trunT\tparticle\tsource\tenergy[MeV]\t";
	for(auto name:nameMap) ofs<<std::to_string(name.first)+"_"+name.second<<"\t"<<massMap[name.first]/g<<"\t";
    ofs<<G4endl;
	ofs.close();
}

RunAction::~RunAction()
{}

G4Run* RunAction::GenerateRun()
{
	// generate run
	fRun = new Run(tetData);
	return fRun;
}


void RunAction::BeginOfRunAction(const G4Run* aRun)
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

	const PrimaryGeneratorAction* primary =
			dynamic_cast<const PrimaryGeneratorAction*>(G4RunManager::GetRunManager()
			->GetUserPrimaryGeneratorAction());
	if(!primary) return;
    if(!primary->GetSourceGenerator()->IsInitialized()){
        G4Exception("Beam Initialization","",FatalErrorInArgument,
                G4String("      Beam direction should be reset after any kind of beam size change").c_str());

    }
	primaryParticle = primary->GetParticleGun()->GetParticleDefinition()->GetParticleName();
	primarySourceName = primary->GetSourceName();
	primaryEnergy = primary->GetParticleGun()->GetParticleEnergy();
    beamArea = primary->GetSourceGenerator()->GetBeamArea();
    fRun->SetPrimary(primaryParticle, primarySourceName, primaryEnergy, beamArea);
}

void RunAction::EndOfRunAction(const G4Run* aRun)
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

	// set doses
	SetDoses();

	// print by G4cout
    PrintResultExternal(G4cout);

	// print by std::ofstream
	std::ofstream ofs(outputFile.c_str(), std::ios::app);
    PrintLineExternal(ofs);
    ofs.close();

	initTimer->Start();
}

void RunAction::SetDoses()
{
	doses.clear();
	EDEPMAP edepMap = *fRun->GetEdepMap();

    for(auto itr : massMap){
		if(itr.first<0) continue;
		G4double meanDose    = edepMap[itr.first].first / numOfEvent;
		G4double squareDoese = edepMap[itr.first].second;
		G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
		G4double relativeE   = sqrt(variance)/meanDose;
        doses[itr.first] = std::make_pair(meanDose/itr.second*1e12, relativeE);
	}
}


void RunAction::PrintResultExternal(std::ostream &out)
{
	// Print run result
	//
	using namespace std;
	EDEPMAP edepMap = *fRun->GetEdepMap();

	out << G4endl
	    << "=======================================================================" << G4endl
	    << " Run #" << runID << " / Number of event processed : "<< numOfEvent     << G4endl
	    << "=======================================================================" << G4endl
		<< " Init time: " << initTimer->GetRealElapsed() << " s / Run time: "<< runTimer->GetRealElapsed()<<" s"<< G4endl
	    << "=======================================================================" << G4endl
	    << setw(27) << "organ ID| "
		<< setw(15) << "Organ Mass (g)"
		<< setw(15) << "Dose (pGy*cm2)"
		<< setw(15) << "Relative Error" << G4endl;

	out.precision(3);

	for(auto itr : massMap){
		if(tetData->DoseWasOrganized()||itr.first<0) out << setw(25) << nameMap[itr.first]<< "| ";
		else                            out << setw(25) << tetData->GetMaterial(itr.first)->GetName()<< "| ";
		out	<< setw(15) << fixed      << itr.second/g;
		out	<< setw(15) << scientific << doses[itr.first].first/(joule/kg)*beamArea/cm2;
		out	<< setw(15) << fixed      << doses[itr.first].second << G4endl;
	}

        out << "=======================================================================" << G4endl << G4endl;
}

void RunAction::PrintLineExternal(std::ostream &out)
{
	// Print run result
	//
	using namespace std;

	out << runID << "\t" <<numOfEvent<<"\t"<< initTimer->GetRealElapsed() << "\t"<< runTimer->GetRealElapsed()<<"\t"
		<< primaryParticle << "\t" <<primarySourceName<< "\t" << primaryEnergy/MeV << "\t";

	for(auto itr:doses){
		out << itr.second.first/(joule/kg) * beamArea/cm2 <<"\t" << itr.second.second << "\t";
    }
    out<<G4endl;
}

std::pair<G4double, G4double> RunAction::PropagateError(std::vector<std::pair<G4double, G4double>> doseVec,
											 std::vector<G4double> ratio)
{
	typedef std::pair<G4double, G4double> VALUE;

	//normalize the ratio
	if(doseVec.size()!=ratio.size()){
		G4Exception("TETRunAction::PropagateError","",JustWarning,
				G4String("      uniform ratio was applied " ).c_str());
		std::vector<G4double> uniform(doseVec.size(),1./(G4int)doseVec.size());
		ratio = uniform;
	}else{
		G4double sum(0.);
		for(size_t i=0;i<ratio.size();i++) sum += ratio[i];
		for(auto &r:ratio) r/=sum;
	}

	//set mean dose
	G4double value(0.);
	for(size_t i=0;i<doseVec.size();i++) value += doseVec[i].first*ratio[i];

	//set error
	G4double error(0.);
	for(size_t i=0;i<doseVec.size();i++){
		if(std::isnan(doseVec[i].second)) continue;
		error += pow(doseVec[i].first*doseVec[i].second*ratio[i],2);
	}
	error = sqrt(error);
	error /= value;

	return VALUE(value, error);
}
