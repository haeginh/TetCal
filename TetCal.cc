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
// External.cc
// \file          Tetcal_v2
// \author        Haegin Han

#include "ActionInitialization.hh"
#include "PhysicsList.hh"
#include "TETDetectorConstruction.hh"
#include "TETModelImport.hh"

#include "G4RunManagerFactory.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"
#include "G4Timer.hh"

void PrintUsage(){
	G4cerr<< "Usage: ./TetCal -m [MACRO] -o [OUTPUT] -p [phantom name] (--usegps)"  <<G4endl;
	G4cerr<< "Example: ./TetCal -m sample.in -o run.out -p ./phantoms/00M" <<G4endl;
}

int main(int argc,char** argv) 
{
	// Read the arguments for batch mode
	//
	G4Timer* initTimer = new G4Timer;
	initTimer->Start();
	G4String macro;
	G4String output("output");
	G4String phantomName;
	G4UIExecutive* ui = 0;
	G4bool useGPS(false);

	for ( G4int i=1; i<argc; i++) {
		// macro file name
		if ( G4String(argv[i]) == "-m" ) {
			macro = argv[++i];
		}
		// output file name
		else if ( G4String(argv[i]) == "-o" ) {
			output = argv[++i];
		}
		// switch for MRCP-AF phantom
		else if ( G4String(argv[i]) == "-p" ) {
			phantomName = argv[++i];
		}
		else if ( G4String(argv[i]) == "--usegps" ) {
			useGPS = true;
		}
		else {
			PrintUsage();
			return 1;
		}
	}

	// print usage when there are more than six arguments
	if (phantomName.empty()){
		G4cout<<"Phantom name is mandatory"<<G4endl;
		PrintUsage();
		return 1;
	}

	// Detect interactive mode (if no macro file name) and define UI session
	//
	G4RunManager* runManager;
	if ( !macro.size() ) {
		ui = new G4UIExecutive(argc, argv);
		runManager = new G4MTRunManager();
	}
	else runManager = G4RunManagerFactory::CreateRunManager();
	// Choose the Random engine
	//
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	G4Random::setTheSeed(time(0));

	// Set a class to import phantom data
	//
	TETModelImport* tetData = new TETModelImport(phantomName, ui);

	// Set mandatory initialisation classes
	//
	// detector construction
	runManager->SetUserInitialization(new TETDetectorConstruction(tetData));
	// physics list
	// runManager->SetUserInitialization(new QBBC);
	runManager->SetUserInitialization(new PhysicsList());
	// user action initialisation
	runManager->SetUserInitialization(new ActionInitialization(tetData, output, initTimer, useGPS));
    
	// Visualization manager
	//
	G4VisManager* visManager = new G4VisExecutive;
	visManager->Initialise();

	// Process macro or start UI session
	//
	G4UImanager* UImanager = G4UImanager::GetUIpointer();
	if ( ! ui ){
		// batch mode
		G4String command = "/control/execute ";
		UImanager->ApplyCommand(command+macro);
	}
	else {
		// interactive mode
		UImanager->ApplyCommand("/control/execute init_vis.mac");
		ui->SessionStart();
		delete ui;
	}

	// Job termination
	//
	delete visManager;
	delete runManager;
}


