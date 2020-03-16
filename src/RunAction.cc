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

#include "RunAction.hh"
#include "G4Timer.hh"
#include <iostream>

RunAction::RunAction(G4String _output, G4Timer* _init)
:fRun(0), numOfEvent(0), runID(0), outputFile(_output), initTimer(_init), runTimer(0),
 primaryDir(1., 0, 0)
{
	if(!isMaster) return;

	runTimer = new G4Timer;
	std::ofstream ofs(outputFile);

	ofs<<"nps\tinitT\trunT\tsource\tdir\t";
	for(G4int i=0;i<2000;i++) ofs<<i<<"\t";
	ofs<<G4endl;
}

RunAction::~RunAction()
{}

G4Run* RunAction::GenerateRun()
{
	// generate run
	fRun = new TETRun();
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

	const TETPrimaryGeneratorAction* primary =
			dynamic_cast<const TETPrimaryGeneratorAction*>(G4RunManager::GetRunManager()
			->GetUserPrimaryGeneratorAction());
	if(!primary) return;
	primarySourceName = primary->GetSourceName();
	primaryDir = primary->GetParticleGun()->GetParticleMomentumDirection();
	fRun->SetPrimary(primarySourceName, primaryDir);

}

void RunAction::EndOfRunAction(const G4Run* aRun)
{
	// print the result only in the Master
	if(!isMaster) return;
	runTimer->Stop();

	// get the run ID
	runID = aRun->GetRunID();

	//get primary info
	primarySourceName  = fRun->GetBeamSourceName();
	primaryDir         = fRun->GetBeamDir();

	// Print the run result by G4cout and std::ofstream
	//

	// print by G4cout
	PrintResult(G4cout);

	// print by std::ofstream
	std::ofstream ofs(outputFile.c_str(), std::ios::app);
	PrintLine(ofs);
	ofs.close();

	initTimer->Start();
}

void RunAction::PrintResult(std::ostream &out)
{
	// Print run result
	//
	using namespace std;
	LENGHBIN lengthBin = fRun->GetLengthBin();

	out << G4endl
	    << "=====================================================================" << G4endl
	    << " Run #" << runID << " / Number of event processed : "<< numOfEvent     << G4endl
	    << "=====================================================================" << G4endl
		<< " Init time: " << initTimer->GetRealElapsed() << " s / Run time: "<< runTimer->GetRealElapsed()<<" s"<< G4endl
	    << "=====================================================================" << G4endl;
	out.precision(3);

	for(auto lb:lengthBin) G4cout<<setw(7)<<lb.first/mm<<" mm   "<<lb.second<<G4endl;

	out << "=====================================================================" << G4endl << G4endl;
}


void RunAction::PrintLine(std::ostream &out)
{
	// Print run result
	//
	using namespace std;
	LENGHBIN lengthBin = fRun->GetLengthBin();

	out << numOfEvent<<"\t"<< initTimer->GetRealElapsed() << "\t"<< runTimer->GetRealElapsed()
		<< "\t" <<primarySourceName<< "\t" << primaryDir << "\t";

	for(G4int i=0;i<2000;i++) out<< lengthBin[i]<<"\t";

	out<<G4endl;
}

