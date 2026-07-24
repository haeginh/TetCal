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
#include "EyeImportanceParallelWorld.hh"
#include "EyeImportanceTuning.hh"
#include "PhysicsList.hh"
#include "ShieldImportanceTuning.hh"
#include "TETDetectorConstruction.hh"
#include "TETModelImport.hh"

#include "G4RunManagerFactory.hh"
#include "G4ParallelWorldPhysics.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#include "G4GeometrySampler.hh"
#include "G4ImportanceBiasing.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4GenericMessenger.hh"

#include "Randomize.hh"
#include "G4Timer.hh"

#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {
std::vector<G4double> ParseNumberList(G4String value, G4double multiplier = 1.);

struct TetCalInputConfig {
	G4String phantomName;
	G4String phantomDataName;
	std::vector<G4String> extraTetNames;
	std::vector<G4String> patientTetNames;
	G4bool eyeImportance = false;
	std::vector<G4double> eyeImportanceRadii = {20.*cm, 10.*cm, 5.*cm, 2.*cm};
	std::vector<G4double> eyeImportanceValues = {1., 2., 5., 10., 20.};
	G4String eyeImportanceParticle = "gamma";
	G4bool eyeImportanceTune = false;
	G4double eyeImportanceTuneTolerance = 0.2;
	G4int eyeImportanceTuneMaxIterations = 10;
	G4bool shieldImportance = false;
	G4double shieldImportanceValue = 1.;
	G4bool shieldImportanceTune = false;
	G4double shieldImportanceTuneTolerance = 0.2;
	G4int shieldImportanceTuneMaxIterations = 10;
	G4String shieldImportanceParticle = "gamma";
	G4int numberOfThreads = 0;

	void SetPhantomName(G4String name) { phantomName = name; }
	void SetPhantomDataName(G4String name) { phantomDataName = name; }
	void AddTetModel(G4String name) { extraTetNames.push_back(name); }
	void AddPatientTetModel(G4String name) { patientTetNames.push_back(name); }
	void ClearTetModels() { extraTetNames.clear(); }
	void ClearPatientTetModels() { patientTetNames.clear(); }
	void SetEyeImportance(G4bool value) { eyeImportance = value; }
	void SetEyeImportanceRadii(G4String value) {
		auto parsed = ParseNumberList(value, cm);
		if(!parsed.empty()) eyeImportanceRadii = parsed;
	}
	void SetEyeImportanceValues(G4String value) {
		auto parsed = ParseNumberList(value);
		if(!parsed.empty()) eyeImportanceValues = parsed;
	}
	void SetEyeImportanceParticle(G4String value) {
		if(!value.empty()) eyeImportanceParticle = value;
	}
	void SetEyeImportanceTune(G4bool value) { eyeImportanceTune = value; }
	void SetEyeImportanceTuneTolerance(G4double value) {
		if(value > 0.) eyeImportanceTuneTolerance = value;
	}
	void SetEyeImportanceTuneMaxIterations(G4int value) {
		if(value > 0) eyeImportanceTuneMaxIterations = value;
	}
	void SetShieldImportance(G4bool value) { shieldImportance = value; }
	void SetShieldImportanceValue(G4double value) {
		if(value > 0.) shieldImportanceValue = value;
	}
	void SetShieldImportanceTune(G4bool value) { shieldImportanceTune = value; }
	void SetShieldImportanceTuneTolerance(G4double value) {
		if(value > 0.) shieldImportanceTuneTolerance = value;
	}
	void SetShieldImportanceTuneMaxIterations(G4int value) {
		if(value > 0) shieldImportanceTuneMaxIterations = value;
	}
	void SetShieldImportanceParticle(G4String value) {
		if(!value.empty()) shieldImportanceParticle = value;
	}
};

G4String Trim(G4String str)
{
	const char* whitespace = " \t\r\n";
	size_t begin = str.find_first_not_of(whitespace);
	if(begin == std::string::npos) return "";
	size_t end = str.find_last_not_of(whitespace);
	return str.substr(begin, end - begin + 1);
}

G4String Unquote(G4String str)
{
	str = Trim(str);
	if(str.size() >= 2 &&
	   ((str.front() == '"' && str.back() == '"') ||
	    (str.front() == '\'' && str.back() == '\''))){
		return str.substr(1, str.size() - 2);
	}
	return str;
}

std::vector<G4double> ParseNumberList(G4String value, G4double multiplier)
{
	std::vector<G4double> result;
	for(auto& c : value){
		if(c == ',') c = ' ';
	}

	std::stringstream ss(value);
	G4double number;
	while(ss >> number){
		result.push_back(number * multiplier);
	}
	return result;
}

G4bool ParseBool(G4String value)
{
	value = Trim(value);
	return value == "1" || value == "true" || value == "True" || value == "TRUE" ||
	       value == "on" || value == "On" || value == "ON";
}

void ReadTetCalInputsFromMacro(G4String macro, TetCalInputConfig& config, G4int depth = 0)
{
	if(macro.empty() || depth > 8) return;

	std::ifstream ifs(macro);
	if(!ifs.is_open()) return;

	G4String line;
	while(std::getline(ifs, line)){
		size_t comment = line.find('#');
		if(comment != std::string::npos) line = line.substr(0, comment);
		line = Trim(line);
		if(line.empty()) continue;

		std::stringstream ss(line);
		G4String command;
		ss >> command;
		G4String value;
		std::getline(ss, value);
		value = Unquote(value);

		if(command == "/tetcal/phantom" || command == "/tetcal/phantomName"){
			config.phantomName = value;
		}
		else if(command == "/tetcal/basePhantom" || command == "/tetcal/basePhantomName"){
			config.phantomDataName = value;
		}
		else if(command == "/tetcal/addTet" || command == "/tetcal/addTetModel"){
			config.extraTetNames.push_back(value);
		}
		else if(command == "/tetcal/addPatientTet" || command == "/tetcal/addPatientTetModel"){
			config.patientTetNames.push_back(value);
		}
		else if(command == "/tetcal/clearTet"){
			config.extraTetNames.clear();
		}
		else if(command == "/tetcal/clearPatientTet"){
			config.patientTetNames.clear();
		}
		else if(command == "/tetcal/eyeImportance"){
			config.eyeImportance = ParseBool(value);
		}
		else if(command == "/tetcal/eyeImportanceRadii"){
			auto parsed = ParseNumberList(value, cm);
			if(!parsed.empty()) config.eyeImportanceRadii = parsed;
		}
		else if(command == "/tetcal/eyeImportanceValues"){
			auto parsed = ParseNumberList(value);
			if(!parsed.empty()) config.eyeImportanceValues = parsed;
		}
		else if(command == "/tetcal/eyeImportanceParticle"){
			if(!value.empty()) config.eyeImportanceParticle = value;
		}
		else if(command == "/tetcal/eyeImportanceTune"){
			config.eyeImportanceTune = ParseBool(value);
		}
		else if(command == "/tetcal/eyeImportanceTuneTolerance"){
			config.eyeImportanceTuneTolerance = std::atof(value.c_str());
		}
		else if(command == "/tetcal/eyeImportanceTuneMaxIterations"){
			config.eyeImportanceTuneMaxIterations = std::atoi(value.c_str());
		}
		else if(command == "/tetcal/shieldImportance"){
			config.shieldImportance = ParseBool(value);
		}
		else if(command == "/tetcal/shieldImportanceValue"){
			config.shieldImportanceValue = std::atof(value.c_str());
		}
		else if(command == "/tetcal/shieldImportanceTune"){
			config.shieldImportanceTune = ParseBool(value);
		}
		else if(command == "/tetcal/shieldImportanceTuneTolerance"){
			config.shieldImportanceTuneTolerance = std::atof(value.c_str());
		}
		else if(command == "/tetcal/shieldImportanceTuneMaxIterations"){
			config.shieldImportanceTuneMaxIterations = std::atoi(value.c_str());
		}
		else if(command == "/tetcal/shieldImportanceParticle"){
			if(!value.empty()) config.shieldImportanceParticle = value;
		}
		else if(command == "/run/numberOfThreads"){
			config.numberOfThreads = std::atoi(value.c_str());
		}
		else if(command == "/control/execute"){
			ReadTetCalInputsFromMacro(value, config, depth + 1);
		}
	}
}

G4bool IsPreInitializationCommand(G4String command)
{
	return command == "/run/numberOfThreads" ||
	       command == "/run/initialize" ||
	       command == "/tetcal/phantom" ||
	       command == "/tetcal/phantomName" ||
	       command == "/tetcal/basePhantom" ||
	       command == "/tetcal/basePhantomName" ||
	       command == "/tetcal/addTet" ||
	       command == "/tetcal/addTetModel" ||
	       command == "/tetcal/addPatientTet" ||
	       command == "/tetcal/addPatientTetModel" ||
	       command == "/tetcal/clearTet" ||
	       command == "/tetcal/clearPatientTet" ||
	       command == "/tetcal/eyeImportance" ||
	       command == "/tetcal/eyeImportanceRadii" ||
	       command == "/tetcal/eyeImportanceValues" ||
	       command == "/tetcal/eyeImportanceParticle" ||
	       command == "/tetcal/eyeImportanceTune" ||
	       command == "/tetcal/eyeImportanceTuneTolerance" ||
	       command == "/tetcal/eyeImportanceTuneMaxIterations" ||
	       command == "/tetcal/shieldImportance" ||
	       command == "/tetcal/shieldImportanceValue" ||
	       command == "/tetcal/shieldImportanceTune" ||
	       command == "/tetcal/shieldImportanceTuneTolerance" ||
	       command == "/tetcal/shieldImportanceTuneMaxIterations" ||
	       command == "/tetcal/shieldImportanceParticle";
}

G4String WritePostInitializationMacro(G4String macro)
{
	if(macro.empty()) return macro;

	std::ifstream ifs(macro);
	if(!ifs.is_open()) return macro;

	G4String filteredMacro = "/private/tmp/tetcal_postinit_" + std::to_string(std::time(nullptr)) + ".mac";
	std::ofstream ofs(filteredMacro);
	if(!ofs.is_open()) return macro;

	G4String line;
	while(std::getline(ifs, line)){
		G4String commandLine = line;
		size_t comment = commandLine.find('#');
		if(comment != std::string::npos) commandLine = commandLine.substr(0, comment);
		commandLine = Trim(commandLine);

		std::stringstream ss(commandLine);
		G4String command;
		ss >> command;

		if(!command.empty() && IsPreInitializationCommand(command)){
			ofs << "# skipped after pre-initialization: " << line << G4endl;
			continue;
		}
		ofs << line << G4endl;
	}

	return filteredMacro;
}
}

void PrintUsage(){
	G4cerr<< "Usage: ./TetCal -m [MACRO] -o [OUTPUT] [-p phantom-node-prefix] [--base-phantom phantom-data-prefix] [--add-tet prefix]... [--add-patient-tet prefix]... [--vis-skin-only] [--eye-importance] (--usegps)"  <<G4endl;
	G4cerr<< "Macro commands: /tetcal/phantom [phantom node prefix], /tetcal/basePhantom [phantom data prefix], /tetcal/addTet [extra TET prefix], /tetcal/addPatientTet [patient TET prefix]" <<G4endl;
	G4cerr<< "                /tetcal/eyeImportance true, /tetcal/eyeImportanceRadii 20 10 5 2, /tetcal/eyeImportanceValues 1 2 5 10 20" <<G4endl;
	G4cerr<< "                /tetcal/eyeImportanceTune true, /tetcal/eyeImportanceTuneTolerance 0.2, /tetcal/eyeImportanceTuneMaxIterations 10" <<G4endl;
	G4cerr<< "                /tetcal/shieldImportance true, /tetcal/shieldImportanceValue 10, /tetcal/shieldImportanceTune true" <<G4endl;
	G4cerr<< "Base phantom reads [phantom node prefix].node; other base files are read from [phantom data prefix].* or MRCP_AM.* in the same directory by default." <<G4endl;
	G4cerr<< "Example: ./TetCal -m sample.in -o run.out" <<G4endl;
}

G4String MakeTunedMacroName(G4String macro)
{
	size_t slash = macro.find_last_of("/\\");
	size_t dot = macro.find_last_of('.');
	if(dot == std::string::npos || (slash != std::string::npos && dot < slash)){
		return macro + "_tuned.in";
	}
	return macro.substr(0, dot) + "_tuned" + macro.substr(dot);
}

G4String MakeTuneDataName(G4String macro)
{
	return macro + ".tuneData";
}

G4bool RunTetCalOnce(int argc,
                    char** argv,
                    G4String macro,
                    G4String output,
                    TetCalInputConfig inputConfig,
                    G4bool useGPS,
                    G4bool visSkinOnly)
{
	G4Timer* initTimer = new G4Timer;
	initTimer->Start();
	G4UIExecutive* ui = 0;

	// Detect interactive mode (if no macro file name) and define UI session
	//
	G4RunManager* runManager;
	if ( !macro.size() ) {
		ui = new G4UIExecutive(argc, argv);
		runManager = new G4MTRunManager();
	}
	else runManager = G4RunManagerFactory::CreateRunManager();

	// Set a class to import phantom data
	//
	TETModelImport* tetData = new TETModelImport(inputConfig.phantomName, inputConfig.phantomDataName, inputConfig.extraTetNames, inputConfig.patientTetNames, ui);

	// Set mandatory initialisation classes
	//
	// detector construction
	TETDetectorConstruction* detector = new TETDetectorConstruction(tetData, useGPS, visSkinOnly);
	EyeImportanceParallelWorld* eyeImportanceWorld = nullptr;
	G4GeometrySampler* eyeImportanceSampler = nullptr;
	static G4int runCounter = 0;
	const G4String eyeImportanceWorldName = "EyeImportanceWorld_" + std::to_string(++runCounter);
	if(inputConfig.eyeImportance){
		eyeImportanceWorld = new EyeImportanceParallelWorld(eyeImportanceWorldName,
		                                                    tetData,
		                                                    inputConfig.eyeImportanceRadii,
		                                                    inputConfig.eyeImportanceValues);
		detector->RegisterParallelWorld(eyeImportanceWorld);
		eyeImportanceSampler = new G4GeometrySampler(eyeImportanceWorldName, inputConfig.eyeImportanceParticle);
		eyeImportanceSampler->SetParallel(true);
	}
	runManager->SetUserInitialization(detector);
	// physics list
	// runManager->SetUserInitialization(new QBBC);
	PhysicsList* physicsList = new PhysicsList();
	if(inputConfig.eyeImportance){
		physicsList->RegisterPhysics(new G4ImportanceBiasing(eyeImportanceSampler, eyeImportanceWorldName));
		physicsList->RegisterPhysics(new G4ParallelWorldPhysics(eyeImportanceWorldName));
	}
	runManager->SetUserInitialization(physicsList);
	// user action initialisation
	runManager->SetUserInitialization(new ActionInitialization(tetData, output, initTimer, useGPS));
    
	// Visualization manager is only needed for interactive sessions.
	G4VisManager* visManager = nullptr;
	if(ui){
		visManager = new G4VisExecutive;
		visManager->Initialise();
	}

	// Process macro or start UI session
	//
	G4UImanager* UImanager = G4UImanager::GetUIpointer();
	G4GenericMessenger tetCalMessenger(&inputConfig, "/tetcal/", "TetCal input control");
	tetCalMessenger.DeclareMethod("phantom", &TetCalInputConfig::SetPhantomName, "Set base phantom node file prefix");
	tetCalMessenger.DeclareMethod("phantomName", &TetCalInputConfig::SetPhantomName, "Set base phantom node file prefix");
	tetCalMessenger.DeclareMethod("basePhantom", &TetCalInputConfig::SetPhantomDataName, "Set base phantom data file prefix for ele/material/dose/RBMnBS/DRF/hands");
	tetCalMessenger.DeclareMethod("basePhantomName", &TetCalInputConfig::SetPhantomDataName, "Set base phantom data file prefix for ele/material/dose/RBMnBS/DRF/hands");
	tetCalMessenger.DeclareMethod("addTet", &TetCalInputConfig::AddTetModel, "Add extra TET file prefix");
	tetCalMessenger.DeclareMethod("addTetModel", &TetCalInputConfig::AddTetModel, "Add extra TET file prefix");
	tetCalMessenger.DeclareMethod("addPatientTet", &TetCalInputConfig::AddPatientTetModel, "Add patient TET file prefix with material IDs remapped to -(100000 + originalID)");
	tetCalMessenger.DeclareMethod("addPatientTetModel", &TetCalInputConfig::AddPatientTetModel, "Add patient TET file prefix with material IDs remapped to -(100000 + originalID)");
	tetCalMessenger.DeclareMethod("clearTet", &TetCalInputConfig::ClearTetModels, "Clear extra TET file prefixes");
	tetCalMessenger.DeclareMethod("clearPatientTet", &TetCalInputConfig::ClearPatientTetModels, "Clear patient TET file prefixes");
	tetCalMessenger.DeclareMethod("eyeImportance", &TetCalInputConfig::SetEyeImportance, "Enable eye-lens centered gamma importance splitting before initialization");
	tetCalMessenger.DeclareMethod("eyeImportanceRadii", &TetCalInputConfig::SetEyeImportanceRadii, "Set eye-importance shell radii in cm");
	tetCalMessenger.DeclareMethod("eyeImportanceValues", &TetCalInputConfig::SetEyeImportanceValues, "Set eye-importance values");
	tetCalMessenger.DeclareMethod("eyeImportanceParticle", &TetCalInputConfig::SetEyeImportanceParticle, "Set particle for eye-importance splitting");
	tetCalMessenger.DeclareMethod("eyeImportanceTune", &TetCalInputConfig::SetEyeImportanceTune, "Iteratively tune eye-importance values in batch mode");
	tetCalMessenger.DeclareMethod("eyeImportanceTuneTolerance", &TetCalInputConfig::SetEyeImportanceTuneTolerance, "Set relative convergence tolerance for eye-importance tuning");
	tetCalMessenger.DeclareMethod("eyeImportanceTuneMaxIterations", &TetCalInputConfig::SetEyeImportanceTuneMaxIterations, "Set maximum eye-importance tuning iterations");
	tetCalMessenger.DeclareMethod("shieldImportance", &TetCalInputConfig::SetShieldImportance, "Enable lead/lead_glass importance splitting");
	tetCalMessenger.DeclareMethod("shieldImportanceValue", &TetCalInputConfig::SetShieldImportanceValue, "Set lead/lead_glass importance value");
	tetCalMessenger.DeclareMethod("shieldImportanceTune", &TetCalInputConfig::SetShieldImportanceTune, "Tune lead/lead_glass importance from enter/exit counts");
	tetCalMessenger.DeclareMethod("shieldImportanceTuneTolerance", &TetCalInputConfig::SetShieldImportanceTuneTolerance, "Set relative convergence tolerance for shield importance tuning");
	tetCalMessenger.DeclareMethod("shieldImportanceTuneMaxIterations", &TetCalInputConfig::SetShieldImportanceTuneMaxIterations, "Set maximum shield importance tuning iterations");
	tetCalMessenger.DeclareMethod("shieldImportanceParticle", &TetCalInputConfig::SetShieldImportanceParticle, "Set particle for shield importance splitting");
	if(inputConfig.numberOfThreads > 0){
		UImanager->ApplyCommand("/run/numberOfThreads " + std::to_string(inputConfig.numberOfThreads));
	}
	if(inputConfig.eyeImportance || inputConfig.shieldImportance){
		runManager->Initialize();
	}
	if(inputConfig.eyeImportance){
		eyeImportanceWorld->CreateImportanceStore();
	}
	if ( ! ui ){
		// batch mode
		G4String command = "/control/execute ";
		G4String macroToExecute = (inputConfig.eyeImportance || inputConfig.shieldImportance) ? WritePostInitializationMacro(macro) : macro;
		UImanager->ApplyCommand(command+macroToExecute);
	}
	else {
		// interactive mode
		UImanager->ApplyCommand("/control/execute init_vis.mac");
		ui->SessionStart();
		delete ui;
	}

	// Job termination
	//
	if(visManager) delete visManager;
	delete runManager;
	delete tetData;
	delete initTimer;
	return true;
}

int main(int argc,char** argv)
{
	// Read the arguments for batch mode
	//
	G4String macro;
	G4String output("output");
	TetCalInputConfig cliConfig;
	G4bool useGPS(false);
	G4bool visSkinOnly(false);

	for ( G4int i=1; i<argc; i++) {
		// macro file name
		if ( G4String(argv[i]) == "-m" ) {
			if(i+1>=argc) {
				PrintUsage();
				return 1;
			}
			macro = argv[++i];
		}
		// output file name
		else if ( G4String(argv[i]) == "-o" ) {
			if(i+1>=argc) {
				PrintUsage();
				return 1;
			}
			output = argv[++i];
		}
		// phantom name
		else if ( G4String(argv[i]) == "-p" ) {
			if(i+1>=argc) {
				PrintUsage();
				return 1;
			}
			cliConfig.phantomName = argv[++i];
		}
		else if ( G4String(argv[i]) == "--add-tet" || G4String(argv[i]) == "-a" ) {
			if(i+1>=argc) {
				PrintUsage();
				return 1;
			}
			cliConfig.extraTetNames.push_back(argv[++i]);
		}
		else if ( G4String(argv[i]) == "--add-patient-tet" ) {
			if(i+1>=argc) {
				PrintUsage();
				return 1;
			}
			cliConfig.patientTetNames.push_back(argv[++i]);
		}
		else if ( G4String(argv[i]) == "--base-phantom" ) {
			if(i+1>=argc) {
				PrintUsage();
				return 1;
			}
			cliConfig.phantomDataName = argv[++i];
		}
		else if ( G4String(argv[i]) == "--usegps" ) {
			useGPS = true;
		}
		else if ( G4String(argv[i]) == "--vis-skin-only" ) {
			visSkinOnly = true;
		}
		else if ( G4String(argv[i]) == "--eye-importance" ) {
			cliConfig.eyeImportance = true;
		}
		else if ( G4String(argv[i]) == "--eye-importance-radii" ) {
			if(i+1>=argc) {
				PrintUsage();
				return 1;
			}
			auto parsed = ParseNumberList(argv[++i], cm);
			if(!parsed.empty()) cliConfig.eyeImportanceRadii = parsed;
		}
		else if ( G4String(argv[i]) == "--eye-importance-values" ) {
			if(i+1>=argc) {
				PrintUsage();
				return 1;
			}
			auto parsed = ParseNumberList(argv[++i]);
			if(!parsed.empty()) cliConfig.eyeImportanceValues = parsed;
		}
		else if ( G4String(argv[i]) == "--eye-importance-particle" ) {
			if(i+1>=argc) {
				PrintUsage();
				return 1;
			}
			cliConfig.eyeImportanceParticle = argv[++i];
		}
		else if ( G4String(argv[i]) == "--eye-importance-tune" ) {
			cliConfig.eyeImportanceTune = true;
		}
		else if ( G4String(argv[i]) == "--shield-importance" ) {
			cliConfig.shieldImportance = true;
		}
		else if ( G4String(argv[i]) == "--shield-importance-value" ) {
			if(i+1>=argc) {
				PrintUsage();
				return 1;
			}
			cliConfig.shieldImportanceValue = std::atof(argv[++i]);
		}
		else if ( G4String(argv[i]) == "--shield-importance-tune" ) {
			cliConfig.shieldImportanceTune = true;
		}
		else {
			PrintUsage();
			return 1;
		}
	}

	TetCalInputConfig initialConfig = cliConfig;
	ReadTetCalInputsFromMacro(macro, initialConfig);

	// print usage when there are more than six arguments
	if (initialConfig.phantomName.empty()){
		G4cout<<"Phantom name is mandatory"<<G4endl;
		PrintUsage();
		return 1;
	}

	// Choose the Random engine
	//
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	G4Random::setTheSeed(time(0));

	if(!macro.empty() && (initialConfig.eyeImportanceTune || initialConfig.shieldImportanceTune)){
		G4String originalMacro = macro;
		G4String tunedMacro = MakeTunedMacroName(originalMacro);
		G4String tuneData = MakeTuneDataName(originalMacro);
		G4String shieldTuneData = originalMacro + ".shieldTuneData";
		G4String currentMacro = originalMacro;
		G4int maxIterations = std::max(initialConfig.eyeImportanceTune ? initialConfig.eyeImportanceTuneMaxIterations : 0,
		                               initialConfig.shieldImportanceTune ? initialConfig.shieldImportanceTuneMaxIterations : 0);

		for(G4int iteration=1; iteration<=maxIterations; iteration++){
			TetCalInputConfig iterationConfig = cliConfig;
			ReadTetCalInputsFromMacro(currentMacro, iterationConfig);
			if(iterationConfig.eyeImportanceTune) iterationConfig.eyeImportance = true;
			if(iterationConfig.shieldImportanceTune) iterationConfig.shieldImportance = true;

			EyeImportanceTuning::Configure(iterationConfig.eyeImportance,
			                               iterationConfig.eyeImportanceRadii,
			                               iterationConfig.eyeImportanceValues,
			                               iterationConfig.eyeImportanceParticle);
			EyeImportanceTuning::ConfigureTuneSession(iterationConfig.eyeImportanceTune,
			                                          currentMacro,
			                                          tunedMacro,
			                                          tuneData,
			                                          iteration,
			                                          iterationConfig.eyeImportanceTuneTolerance);
			ShieldImportanceTuning::Configure(iterationConfig.shieldImportance,
			                                  iterationConfig.shieldImportanceValue,
			                                  iterationConfig.shieldImportanceParticle);
			ShieldImportanceTuning::ConfigureTuneSession(iterationConfig.shieldImportanceTune,
			                                             currentMacro,
			                                             tunedMacro,
			                                             shieldTuneData,
			                                             iteration,
			                                             iterationConfig.shieldImportanceTuneTolerance);

			G4cout << "Starting importance tuning iteration " << iteration
			       << " using macro " << currentMacro << G4endl;
			RunTetCalOnce(argc, argv, currentMacro, output, iterationConfig, useGPS, visSkinOnly);

			G4bool eyeDone = !iterationConfig.eyeImportanceTune || EyeImportanceTuning::IsTuneConverged();
			G4bool shieldDone = !iterationConfig.shieldImportanceTune || ShieldImportanceTuning::IsTuneConverged();
			if(eyeDone && shieldDone){
				G4cout << "Importance tuning converged after " << iteration
				       << " iteration(s). Tuned macro: " << tunedMacro
				       << ", eye tune data: " << tuneData
				       << ", shield tune data: " << shieldTuneData << G4endl;
				break;
			}

			currentMacro = tunedMacro;
			if(iteration == maxIterations){
				G4cout << "Importance tuning reached the maximum iteration count ("
				       << maxIterations << "). Last tuned macro: " << tunedMacro
				       << ", eye tune data: " << tuneData
				       << ", shield tune data: " << shieldTuneData << G4endl;
			}
		}
		return 0;
	}

	EyeImportanceTuning::Configure(initialConfig.eyeImportance,
	                               initialConfig.eyeImportanceRadii,
	                               initialConfig.eyeImportanceValues,
	                               initialConfig.eyeImportanceParticle);
	EyeImportanceTuning::ConfigureTuneSession(false, "", "", "", 0, initialConfig.eyeImportanceTuneTolerance);
	ShieldImportanceTuning::Configure(initialConfig.shieldImportance,
	                                  initialConfig.shieldImportanceValue,
	                                  initialConfig.shieldImportanceParticle);
	ShieldImportanceTuning::ConfigureTuneSession(false, "", "", "", 0, initialConfig.shieldImportanceTuneTolerance);
	RunTetCalOnce(argc, argv, macro, output, initialConfig, useGPS, visSkinOnly);
	return 0;
}
