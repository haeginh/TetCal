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
// \update
// \  - add internal source
// \  -
// \  -
// \  -
// \  -

#include "TETDetectorConstruction.hh"
#include "TETModelImport.hh"
#include "TETPhysicsList.hh"
#include "TETActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#ifdef G4UI_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4VIS_USE
#include "G4UIExecutive.hh"
#endif

#include "Randomize.hh"
#include "G4Timer.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"

void PrintUsage(){
	G4cerr<< "Usage: ./TetCal -m [MACRO] -o [OUTPUT] -p [phantom name]"  <<G4endl;
	G4cerr<< "Example: ./TetCal -m sample.in -o run.out -p ./phantoms/00M" <<G4endl;
}

int main(int argc,char** argv) 
{
	// Read the arguments for batch mode
	//
	G4Timer* initTimer = new G4Timer;
	initTimer->Start();
	G4String macro;
	G4String output;
	G4String phantomName;
	G4UIExecutive* ui = 0;

	for ( G4int i=1; i<argc; i++ ) {
		// macro file name
		if ( G4String(argv[i]) == "-m" ) {
			macro = argv[i+1];
			i++;
		}
		// output file name
		else if ( G4String(argv[i]) == "-o" ) {
			output = argv[i+1];
			i++;
		}
		// switch for MRCP-AF phantom
		else if ( G4String(argv[i]) == "-p" ) {
			phantomName = argv[i+1];
			i++;
		}
		else {
			PrintUsage();
			return 1;
		}
	}

	// print usage when there are more than six arguments
	if ( argc>7 || macro.empty() || phantomName.empty()){
		PrintUsage();
		return 1;
	}

	// Detect interactive mode (if no macro file name) and define UI session
	//
	if ( !macro.size() ) {
#ifdef G4UI_USE
		ui = new G4UIExecutive(argc, argv, "csh");
#else
		G4cerr<<"ERROR: Interactive mode is not available. Please provide macro file."<<G4endl;
		return 1;
#endif
	}
	// default output file name
	else if ( !output.size() ) output = macro + ".out";

	// Choose the Random engine
	//
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	G4Random::setTheSeed(time(0));

	// Construct the default run manager
	//
#ifdef G4MULTITHREADED
	G4MTRunManager * runManager = new G4MTRunManager;
	// set the default number of threads as one
	runManager->SetNumberOfThreads(1);
#else
	G4RunManager * runManager = new G4RunManager;
#endif

	// Set a class to import phantom data
	//
	TETModelImport* tetData = new TETModelImport(phantomName, ui);

	// Set mandatory initialisation classes
	//
	// detector construction
	runManager->SetUserInitialization(new TETDetectorConstruction(tetData));
	// physics list
	G4PhysListFactory factory;
//	G4VModularPhysicsList* physList = factory.GetReferencePhysList("QGSP_BIC_LIV");
//	runManager->SetUserInitialization(physList);
	runManager->SetUserInitialization(new TETPhysicsList());
	// user action initialisation
	runManager->SetUserInitialization(new TETActionInitialization(tetData, output, initTimer));
    
#ifdef G4VIS_USE
	// Visualization manager
	//
	G4VisManager* visManager = new G4VisExecutive;
	visManager->Initialise();
#endif

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
#ifdef G4VIS_USE
		UImanager->ApplyCommand("/control/execute init_vis.mac");
#endif
		ui->SessionStart();
#ifdef G4VIS_USE
		delete visManager;
#endif
		delete ui;
	}

	// Job termination
	//
	delete runManager;
}


