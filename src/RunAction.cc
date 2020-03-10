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
#include "../include/RunAction.hh"

RunAction::RunAction(TETModelImport* _tetData, G4String _output, G4Timer* _init)
:tetData(_tetData), fRun(0), numOfEvent(0), runID(0), outputFile(_output), initTimer(_init), runTimer(0),
 primaryEnergy(-1.), beamArea(-1.), isExternal(true)
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

//	nameMap[-6] = "RBM(DRF)_test"; nameMap[-5] = "BS(DRF)_test";
	nameMap[-4] = "RBM(DRF)"; nameMap[-3] = "BS(DRF)";
	nameMap[-2] = "RBM"     ; nameMap[-1] = "BS"     ;

	ofs<<"[External: pGycm2 / Internal: SAF (kg-1)]"<<G4endl;
	ofs<<"run#\tnps\tinitT\trunT\tparticle\tsource\tenergy[MeV]\t";
	for(auto name:nameMap) ofs<<std::to_string(name.first)+"_"+name.second<<"\t"<<massMap[name.first]/g<<"\t";
	ofs<<"eff. dose (DRF)"<<"\t\t"<< "eff. dose"<<G4endl;
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
	primaryParticle = primary->GetParticleGun()->GetParticleDefinition()->GetParticleName();
	primarySourceName = primary->GetSourceName();
	primaryEnergy = primary->GetParticleGun()->GetParticleEnergy();
	beamArea = primary->GetExternalBeamGenerator()->GetBeamArea();
	isExternal = primary-> GetSourceGenerator()->IsExternal();
	fRun->SetPrimary(primaryParticle, primarySourceName, primaryEnergy, beamArea, isExternal);

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
	isExternal      = fRun->GetIsExternal();

	G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle(primaryParticle);
	weight = GetRadiationWeighting(particle, primaryEnergy);
	// Print the run result by G4cout and std::ofstream
	//

	// set doses
	SetDoses();
	if(isExternal) SetEffectiveDose();

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

void RunAction::SetDoses()
{
	doses.clear();
	EDEPMAP edepMap = *fRun->GetEdepMap();
	for(G4int i=-4;i<0;i++){
		G4double meanDose = edepMap[i].first / numOfEvent;
		G4double squareDoese = edepMap[i].second;
		G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
		G4double relativeE   = sqrt(variance)/meanDose;
		doses[i] = std::make_pair(meanDose, relativeE);
	}
	if(isExternal){
		doses[-1].first *= 1e12;
		doses[-2].first *= 1e12;
	}
	for(auto itr : massMap){
		if(itr.first<0) continue;
		G4double meanDose    = edepMap[itr.first].first  / itr.second / numOfEvent;
		G4double squareDoese = edepMap[itr.first].second / (itr.second*itr.second);
		G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
		G4double relativeE   = sqrt(variance)/meanDose;
		if(isExternal) doses[itr.first] = std::make_pair(meanDose*1e12, relativeE);
		else           doses[itr.first] = std::make_pair(meanDose, relativeE);
	}

	for(auto itr:doses) itr.second.first *= weight;
}

void RunAction::SetEffectiveDose()
{
	G4int group1_RBM = -2;
	G4int group1_RBM_DRF = -4;
	G4int group4_BS = -1;
	G4int group4_BS_DRF = -3;
	G4int ET1(21), ET2(22);
	std::vector<G4int> group1 = {2, 4, 5, 7};
	std::vector<G4int> group2 = {8};
	std::vector<G4int> group3 = {9, 11, 13, 14};
	std::vector<G4int> group4 = {16, 17, 19};
	std::vector<G4int> remainder = {20, 25, 26, 27, 28, 29, 30, 31, 32, 33, 35, 36};

	std::vector<std::pair<G4double, G4double>> effDoseComp;
	effDoseComp.push_back(doses[group1_RBM]);
	for(G4int idx:group1)	effDoseComp.push_back(doses[idx]);
	for(G4int idx:group2)	effDoseComp.push_back(doses[idx]);
	for(G4int idx:group3)	effDoseComp.push_back(doses[idx]);
	for(G4int idx:group4)	effDoseComp.push_back(doses[idx]);
	effDoseComp.push_back(doses[group4_BS]);
	for(G4int idx:remainder)	effDoseComp.push_back(doses[idx]);
	effDoseComp.push_back(doses[ET1]); effDoseComp.push_back(doses[ET2]);

	std::vector<G4double> ratios;
    for(size_t i=0;i<group1.size()+1;i++) ratios.push_back(0.12);
    for(size_t i=0;i<group2.size();i++) ratios.push_back(0.08);
    for(size_t i=0;i<group3.size();i++) ratios.push_back(0.04);
    for(size_t i=0;i<group4.size()+1;i++) ratios.push_back(0.01);
    for(size_t i=0;i<remainder.size();i++) ratios.push_back(0.12/(G4double)(remainder.size()+1));
    ratios.push_back(0.12/(G4double)(remainder.size()+1)*0.001);
    ratios.push_back(0.12/(G4double)(remainder.size()+1)*0.999);

    std::vector<std::pair<G4double, G4double>> effDoseComp_DRF(effDoseComp);
	effDoseComp_DRF[0] = doses[group1_RBM_DRF];
	effDoseComp_DRF[13] = doses[group4_BS_DRF];

	G4double sum(0.);
	for(auto r:ratios) sum += r;
	G4cout<<"Sum of the ratios -->"<<sum<<G4endl;

	effective = PropagateError(effDoseComp, ratios);
	effective_DRF = PropagateError(effDoseComp_DRF, ratios);
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

//	for(G4int i=-4;i<0;i++){
//		out << setw(25) << nameMap[i] << "| ";
//		out << setw(30) << scientific << doses[i].first/(joule/kg)*beamArea/cm2<< setw(15) << fixed << doses[i].second << G4endl;
//	}

	for(auto itr : massMap){
		if(tetData->DoseWasOrganized()||itr.first<0) out << setw(25) << nameMap[itr.first]<< "| ";
		else                            out << setw(25) << tetData->GetMaterial(itr.first)->GetName()<< "| ";
		out	<< setw(15) << fixed      << itr.second/g;
		out	<< setw(15) << scientific << doses[itr.first].first/(joule/kg)*beamArea/cm2;
		out	<< setw(15) << fixed      << doses[itr.first].second << G4endl;
	}

	//effective dose
	out << setw(25) << "eff. dose (DRF)" << "| ";
	out	<< setw(15) << " "                ;
	out	<< setw(15) << scientific << effective_DRF.first/(joule/kg)*beamArea/cm2;
	out	<< setw(15) << fixed      << effective_DRF.second << G4endl;

	out << setw(25) << "eff. dose" << "| ";
	out	<< setw(15) << " "                ;
	out	<< setw(15) << scientific << effective.first/(joule/kg)*beamArea/cm2;
	out	<< setw(15) << fixed      << effective.second << G4endl;

	out << "=======================================================================" << G4endl << G4endl;
}

void RunAction::PrintResultInternal(std::ostream &out)
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
		<< setw(15) << "SAF (kg-1)"
		<< setw(15) << "Relative Error" << G4endl;

	out.precision(3);

	for(auto itr : massMap){
		if(tetData->DoseWasOrganized()||itr.first<0) out << setw(25) << nameMap[itr.first]<< "| ";
		else                            out << setw(25) << tetData->GetMaterial(itr.first)->GetName()<< "| ";
		out	<< setw(15) << fixed      << itr.second/g;
		out	<< setw(15) << scientific << doses[itr.first].first/primaryEnergy/(1./kg);
		out	<< setw(15) << fixed      << doses[itr.first].second << G4endl;
	}
	out << "=======================================================================" << G4endl << G4endl;
}

void RunAction::PrintLineExternal(std::ostream &out)
{
	// Print run result
	//
	using namespace std;
	EDEPMAP edepMap = *fRun->GetEdepMap();

	out << runID << "\t" <<numOfEvent<<"\t"<< initTimer->GetRealElapsed() << "\t"<< runTimer->GetRealElapsed()<<"\t"
		<< primaryParticle << "\t" <<primarySourceName<< "\t" << primaryEnergy/MeV << "\t";

	for(auto itr:doses){
		out << itr.second.first/(joule/kg) * beamArea/cm2 <<"\t" << itr.second.second << "\t";
	}
	out<<effective_DRF.first/(joule/kg) * beamArea/cm2<< "\t" <<effective_DRF.second <<"\t";
	out<<effective.first/(joule/kg) * beamArea/cm2<< "\t" <<effective.second ;
	out<<G4endl;
}

void RunAction::PrintLineInternal(std::ostream &out)
{
	// Print run result
	//
	using namespace std;
	EDEPMAP edepMap = *fRun->GetEdepMap();

	out << runID << "\t" <<numOfEvent<<"\t"<< initTimer->GetRealElapsed() << "\t"<< runTimer->GetRealElapsed()<<"\t"
		<< primaryParticle << "\t" <<primarySourceName<< "\t" << primaryEnergy/MeV << "\t";

	for(auto itr:doses){
		out << itr.second.first/primaryEnergy/(1./kg) <<"\t" << itr.second.second << "\t";
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

G4double RunAction::GetRadiationWeighting(G4ParticleDefinition* _particle, G4double _energy)
{
	G4double weightingFactor = 1.0; // for Gamma and Electron, Note that Muons need to be considered later.

	if(_particle == G4Proton::Proton()) { //charged pions need to be considered later.
		weightingFactor = 2.0;
	}
	else if(_particle == G4Alpha::Alpha()) { // fission fragments and heavy ions need to be considered later.
		weightingFactor = 20.0;
	}
	else if(_particle == G4Neutron::Neutron()) { //neutron
		if( _energy < 1.0) {// under 1 MeV
			weightingFactor = 2.5 + 18.2*exp(-(pow((log(_energy)),2))/6);
		}
		else if(_energy <= 50) { // From 1 MeV to 50 MeV
			weightingFactor = 5.0 + 17.0*exp(-(pow((log(2.0*_energy)),2))/6);
		}
		else {// more than 50 MeV
			weightingFactor = 2.5 + 3.25*exp(-(pow((log(0.04*_energy)),2))/6);
		}
	}

	return weightingFactor;

}


