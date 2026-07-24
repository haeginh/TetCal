#include "ShieldImportanceTuning.hh"

#include "G4ios.hh"
#include "G4AutoLock.hh"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>

G4bool ShieldImportanceTuning::enabled = false;
G4bool ShieldImportanceTuning::tuneEnabled = false;
G4bool ShieldImportanceTuning::tuneConverged = false;
G4double ShieldImportanceTuning::importance = 1.;
G4String ShieldImportanceTuning::particleName = "gamma";
G4String ShieldImportanceTuning::tuneInputMacro;
G4String ShieldImportanceTuning::tuneOutputMacro;
G4String ShieldImportanceTuning::tuneDataFile;
G4int ShieldImportanceTuning::tuneIteration = 0;
G4double ShieldImportanceTuning::tuneTolerance = 0.2;
G4double ShieldImportanceTuning::enter = 0.;
G4double ShieldImportanceTuning::weightedEnter = 0.;
G4double ShieldImportanceTuning::exit = 0.;
G4double ShieldImportanceTuning::weightedExit = 0.;

namespace {
G4Mutex shieldImportanceMutex = G4MUTEX_INITIALIZER;
}

namespace {
G4String TrimLine(G4String str)
{
	const char* whitespace = " \t\r\n";
	size_t begin = str.find_first_not_of(whitespace);
	if(begin == std::string::npos) return "";
	size_t end = str.find_last_not_of(whitespace);
	return str.substr(begin, end - begin + 1);
}

G4String FirstCommand(G4String line)
{
	size_t comment = line.find('#');
	if(comment != std::string::npos) line = line.substr(0, comment);
	line = TrimLine(line);
	std::stringstream ss(line);
	G4String command;
	ss >> command;
	return command;
}

G4String ImportanceLine(G4double value)
{
	std::ostringstream os;
	os << "/tetcal/shieldImportanceValue " << std::setprecision(8) << value;
	return os.str();
}

G4double RelativeDifference(G4double current, G4double suggested)
{
	if(current == 0.) return suggested == 0. ? 0. : 1.e30;
	return std::fabs(suggested - current) / std::fabs(current);
}
}

void ShieldImportanceTuning::Configure(G4bool _enabled,
                                       G4double _importance,
                                       const G4String& _particleName)
{
	enabled = _enabled;
	importance = std::max(1., _importance);
	particleName = _particleName.empty() ? "gamma" : _particleName;
}

void ShieldImportanceTuning::ConfigureTuneSession(G4bool _enabled,
                                                  const G4String& inputMacro,
                                                  const G4String& tunedMacro,
                                                  const G4String& tuneData,
                                                  G4int iteration,
                                                  G4double tolerance)
{
	tuneEnabled = _enabled;
	tuneConverged = false;
	tuneInputMacro = inputMacro;
	tuneOutputMacro = tunedMacro;
	tuneDataFile = tuneData;
	tuneIteration = iteration;
	tuneTolerance = tolerance;
}

G4bool ShieldImportanceTuning::IsEnabled() { return enabled; }
G4bool ShieldImportanceTuning::IsTuneEnabled() { return tuneEnabled; }
G4bool ShieldImportanceTuning::IsTuneConverged() { return tuneConverged; }
G4double ShieldImportanceTuning::GetImportance() { return importance; }
const G4String& ShieldImportanceTuning::GetParticleName() { return particleName; }

void ShieldImportanceTuning::ResetCounts()
{
	G4AutoLock lock(&shieldImportanceMutex);
	enter = 0.;
	weightedEnter = 0.;
	exit = 0.;
	weightedExit = 0.;
}

void ShieldImportanceTuning::AddEnter(G4double weight)
{
	if(!enabled) return;
	G4AutoLock lock(&shieldImportanceMutex);
	enter += 1.;
	weightedEnter += weight;
}

void ShieldImportanceTuning::AddExit(G4double weight)
{
	if(!enabled) return;
	G4AutoLock lock(&shieldImportanceMutex);
	exit += 1.;
	weightedExit += weight;
}

G4double ShieldImportanceTuning::GetEnter()
{
	G4AutoLock lock(&shieldImportanceMutex);
	return enter;
}

G4double ShieldImportanceTuning::GetWeightedEnter()
{
	G4AutoLock lock(&shieldImportanceMutex);
	return weightedEnter;
}

G4double ShieldImportanceTuning::GetExit()
{
	G4AutoLock lock(&shieldImportanceMutex);
	return exit;
}

G4double ShieldImportanceTuning::GetWeightedExit()
{
	G4AutoLock lock(&shieldImportanceMutex);
	return weightedExit;
}

G4double ShieldImportanceTuning::RecommendImportance(G4double weightedEnter, G4double weightedExit)
{
	if(weightedEnter <= 0. || weightedExit <= 0.) return importance;
	return std::max(1., weightedEnter / weightedExit);
}

void ShieldImportanceTuning::PrintSummary(std::ostream& out,
                                          G4double enter,
                                          G4double weightedEnter,
                                          G4double exit,
                                          G4double weightedExit)
{
	if(!enabled) return;
	G4double suggested = RecommendImportance(weightedEnter, weightedExit);
	out << G4endl
	    << "=======================================================================" << G4endl
	    << " Shield importance tuning summary for particle: " << particleName << G4endl
	    << " Materials: lead(-1), lead_glass(-2), syringe_shield_lead(-7)" << G4endl
	    << " TrackEnter | WeightedEnter | TrackExit | WeightedExit | CurrentImp | SuggestedImp" << G4endl;
	out.setf(std::ios::fixed, std::ios::floatfield);
	out.precision(3);
	out << std::setw(10) << enter
	    << " | " << std::setw(13) << weightedEnter
	    << " | " << std::setw(9) << exit
	    << " | " << std::setw(12) << weightedExit
	    << " | " << std::setw(10) << importance
	    << " | " << std::setw(12) << suggested
	    << G4endl
	    << " Recommended macro:" << G4endl
	    << ImportanceLine(suggested) << G4endl
	    << "=======================================================================" << G4endl
	    << G4endl;
}

void ShieldImportanceTuning::WriteTunedMacroAndData(G4double enter,
                                                    G4double weightedEnter,
                                                    G4double exit,
                                                    G4double weightedExit)
{
	if(!enabled || !tuneEnabled || tuneInputMacro.empty() || tuneOutputMacro.empty()) return;

	G4double suggested = RecommendImportance(weightedEnter, weightedExit);
	G4double relative = RelativeDifference(importance, suggested);
	tuneConverged = relative <= tuneTolerance;

	std::ifstream ifs(tuneOutputMacro);
	if(!ifs.is_open()){
		ifs.open(tuneInputMacro);
	}
	if(ifs.is_open()){
		G4String tmpTarget = tuneOutputMacro + ".shield.tmp";
		std::ofstream ofs(tmpTarget);
		G4bool replaced = false;
		G4bool inserted = false;
		G4String line;
		while(std::getline(ifs, line)){
			G4String command = FirstCommand(line);
			if(command == "/tetcal/shieldImportanceValue"){
				ofs << ImportanceLine(suggested) << G4endl;
				replaced = true;
				continue;
			}
			ofs << line << G4endl;
			if(!replaced && !inserted && command == "/tetcal/shieldImportance"){
				ofs << ImportanceLine(suggested) << G4endl;
				inserted = true;
			}
		}
		if(!replaced && !inserted){
			ofs << "/tetcal/shieldImportance true" << G4endl;
			ofs << ImportanceLine(suggested) << G4endl;
		}
		ofs.close();
		if(ifs.good() || ifs.eof()){
			std::remove(tuneOutputMacro.c_str());
			std::rename(tmpTarget.c_str(), tuneOutputMacro.c_str());
		}
	}

	if(!tuneDataFile.empty()){
		std::ofstream data;
		if(tuneIteration <= 1){
			data.open(tuneDataFile);
			data << "iteration\ttrack_enter\tweighted_enter\ttrack_exit\tweighted_exit\tcurrent_importance\tsuggested_importance\trelative_difference\tconverged" << G4endl;
		}
		else{
			data.open(tuneDataFile, std::ios::app);
		}
		data << tuneIteration << "\t"
		     << enter << "\t"
		     << weightedEnter << "\t"
		     << exit << "\t"
		     << weightedExit << "\t"
		     << importance << "\t"
		     << suggested << "\t"
		     << relative << "\t"
		     << (tuneConverged ? "yes" : "no") << G4endl;
	}

	G4cout << "Shield importance tune iteration " << tuneIteration
	       << ": wrote " << tuneOutputMacro
	       << ", relative change = " << relative
	       << ", converged = " << (tuneConverged ? "yes" : "no")
	       << G4endl;
}
