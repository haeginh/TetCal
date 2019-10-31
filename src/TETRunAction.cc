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
	std::ofstream ofs(outputFile);

	if(tetData->DoseWasOrganized()) massMap = tetData->GetDoseMassMap();
	else	                        massMap = tetData->GetMassMap();

	ofs<<"[External: pGycm2 / Internal: SAF (kg-1)]"<<G4endl;
	ofs<<"run#\tnps\tinitT\trunT\tparticle\tsource\tenergy[MeV]\tRBM_homo\t\tBS_homo\t\tRBM_DRF\t\tBS_DRF\t\t";
	for(auto itr : massMap)
		if(tetData->DoseWasOrganized())	ofs<<std::to_string(itr.first)+"_"+tetData->GetDoseName(itr.first)<<"\t"<<itr.second/g<<"\t";
		else                            ofs<<std::to_string(itr.first)+"_"+tetData->GetMaterial(itr.first)->GetName()<<"\t"<<itr.second/g<<"\t";
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

void TETRunAction::SetDoses()
{
	doseValues.clear(); doseErrors.clear();
	EDEPMAP edepMap = *fRun->GetEdepMap();
	for(G4int i=-4;i<0;i++){
		G4double meanDose = edepMap[i].first / numOfEvent;
		G4double squareDoese = edepMap[i].second;
		G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
		G4double relativeE   = sqrt(variance)/meanDose;
		doseValues.push_back(meanDose);
		doseErrors.push_back(relativeE);
	}

	for(auto itr : massMap){
		G4double meanDose    = edepMap[itr.first].first  / itr.second / numOfEvent;
		G4double squareDoese = edepMap[itr.first].second / (itr.second*itr.second);
		G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
		G4double relativeE   = sqrt(variance)/meanDose;
		doseValues.push_back(meanDose);
		doseErrors.push_back(relativeE);
	}
}

void TETRunAction::SetEffectiveDose()
{

}

void TETRunAction::PrintResultExternal(std::ostream &out)
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
	for(i=0;i<4;i++){
		if(i==0) out << setw(27) << "RBM_homo| ";
		if(i==1) out << setw(27) << "BS_homo| ";
		if(i==2) out << setw(27) << "RBM_DRF| ";
		if(i==3) out << setw(27) << "BS_DRF| ";
		out << setw(30) << scientific << doseValues[i]/(joule/kg)*beamArea/cm2<< setw(15) << fixed << doseErrors[i] << G4endl;
	}

	for(auto itr : massMap){
		if(tetData->DoseWasOrganized()) out << setw(25) << tetData->GetDoseName(itr.first)<< "| ";
		else                            out << setw(25) << tetData->GetMaterial(itr.first)->GetName()<< "| ";
		out	<< setw(15) << fixed      << itr.second/g;
		out	<< setw(15) << scientific << doseValues[i]/(joule/kg)*beamArea/cm2;
		out	<< setw(15) << fixed      << doseErrors[i++] << G4endl;
	}

	out << "=====================================================================" << G4endl << G4endl;
}

void TETRunAction::PrintResultInternal(std::ostream &out)
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
	for(i=0;i<4;i++){
		if(i==0) out << setw(27) << "RBM_homo| ";
		if(i==1) out << setw(27) << "BS_homo| ";
		if(i==2) out << setw(27) << "RBM_DRF| ";
		if(i==3) out << setw(27) << "BS_DRF| ";
		out << setw(30) << scientific << doseValues[i]/primaryEnergy/(1./kg)<< setw(15) << fixed << doseErrors[i] << G4endl;
	}

	for(auto itr : massMap){
		if(tetData->DoseWasOrganized()) out << setw(25) << tetData->GetDoseName(itr.first)<< "| ";
		else                            out << setw(25) << tetData->GetMaterial(itr.first)->GetName()<< "| ";
		out	<< setw(15) << fixed      << itr.second/g;
		out	<< setw(15) << scientific << doseValues[i]/primaryEnergy/(1./kg);
		out	<< setw(15) << fixed      << doseErrors[i++] << G4endl;
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

	for(size_t i=0;i<doseValues.size();i++){
		out << doseValues[i]*1e12/(joule/kg) * beamArea/cm2 <<"\t" << doseErrors[i] << "\t";
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

	for(size_t i=0;i<doseValues.size();i++){
		out << doseValues[i]/primaryEnergy/(1./kg) <<"\t" << doseErrors[i] << "\t";
	}
	out<<G4endl;
}

std::pair<G4double, G4double> PropagateError(std::vector<G4double> doses,
										     std::vector<G4double> errors,
											 std::vector<G4double> ratio)
{
	typedef std::pair<G4double, G4double> VALUE;

	if(doses.size()!=errors.size()){
		G4Exception("TETRunAction::PropagateError","",JustWarning,
				G4String("      wrong input " ).c_str());
		return VALUE(0, 0);
	}

	//normalize the ratio
	if(doses.size()!=ratio.size()){
		G4Exception("TETRunAction::PropagateError","",JustWarning,
				G4String("      uniform ratio was applied " ).c_str());
		std::vector<G4double> uniform(doses.size(),1./(G4int)doses.size());
		ratio = uniform;
	}else{
		G4double sum(0.);
		for(size_t i=0;i<ratio.size();i++) sum += ratio[i];
		for(auto &r:ratio) r/=sum;
	}

	//set mean dose
	G4double value(0.);
	for(size_t i=0;i<doses.size();i++) value += doses[i]*ratio[i];

	//set error
	G4double error(0.);
	for(size_t i=0;i<doses.size();i++) error += pow(doses[i]*errors[i]*ratio[i],2);
	error = sqrt(error);

	return VALUE(value, error);
}
