#include "EyeImportanceTuning.hh"

#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>

G4bool EyeImportanceTuning::enabled = false;
G4bool EyeImportanceTuning::tuneEnabled = false;
G4bool EyeImportanceTuning::tuneConverged = false;
std::vector<G4double> EyeImportanceTuning::radii;
std::vector<G4double> EyeImportanceTuning::importances;
G4String EyeImportanceTuning::particleName = "gamma";
G4String EyeImportanceTuning::tuneInputMacro;
G4String EyeImportanceTuning::tuneOutputMacro;
G4String EyeImportanceTuning::tuneDataFile;
G4int EyeImportanceTuning::tuneIteration = 0;
G4double EyeImportanceTuning::tuneTolerance = 0.2;

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

G4String ImportanceLine(const std::vector<G4double>& values)
{
	std::ostringstream os;
	os << "/tetcal/eyeImportanceValues ";
	for(auto value : values) os << std::setprecision(8) << value << " ";
	return os.str();
}

G4double RelativeDifference(G4double current, G4double suggested)
{
	if(current == 0.) return suggested == 0. ? 0. : 1.e30;
	return std::fabs(suggested - current) / std::fabs(current);
}
}

void EyeImportanceTuning::Configure(G4bool _enabled,
                                    const std::vector<G4double>& _radii,
                                    const std::vector<G4double>& _importances,
                                    const G4String& _particleName)
{
	enabled = _enabled;
	radii = _radii;
	importances = _importances;
	particleName = _particleName;
	std::sort(radii.begin(), radii.end(), std::greater<G4double>());
}

void EyeImportanceTuning::ConfigureTuneSession(G4bool _enabled,
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

G4bool EyeImportanceTuning::IsEnabled()
{
	return enabled;
}

G4bool EyeImportanceTuning::IsTuneEnabled()
{
	return tuneEnabled;
}

G4bool EyeImportanceTuning::IsTuneConverged()
{
	return tuneConverged;
}

const std::vector<G4double>& EyeImportanceTuning::GetRadii()
{
	return radii;
}

const std::vector<G4double>& EyeImportanceTuning::GetImportances()
{
	return importances;
}

const G4String& EyeImportanceTuning::GetParticleName()
{
	return particleName;
}

const G4String& EyeImportanceTuning::GetTunedMacroName()
{
	return tuneOutputMacro;
}

std::vector<G4double> EyeImportanceTuning::RecommendImportances(const std::map<G4int, G4double>& weightedCounts,
                                                                G4int nEvents,
                                                                G4double maxRatio)
{
	std::vector<G4double> recommended;
	recommended.push_back(1.);

	if(nEvents <= 0) return recommended;
	if(maxRatio < 1.) maxRatio = 1.;

	for(G4int i=0; i<(G4int)radii.size(); i++){
		G4int cellID = i + 1;
		G4double raw = recommended.back();
		auto count = weightedCounts.find(cellID);
		if(count != weightedCounts.end() && count->second > 0.){
			raw = (G4double)nEvents / count->second;
		}

		G4double lower = recommended.back();
		G4double upper = recommended.back() * maxRatio;
		recommended.push_back(std::min(std::max(raw, lower), upper));
	}

	return recommended;
}

void EyeImportanceTuning::PrintSummary(std::ostream& out,
                                       const std::map<G4int, G4double>& counts,
                                       const std::map<G4int, G4double>& weightedCounts,
                                       G4int nEvents)
{
	if(!enabled || radii.empty()) return;

	auto recommended = RecommendImportances(weightedCounts, nEvents);

	out << G4endl
	       << "=======================================================================" << G4endl
	       << " Eye importance tuning summary for particle: " << particleName << G4endl
	       << " Target track population per shell: " << nEvents << G4endl
	       << " Cell | Radius(cm) | TrackEnter | WeightedEnter | CurrentImp | SuggestedImp" << G4endl;

	out.setf(std::ios::fixed, std::ios::floatfield);
	out.precision(3);

	for(G4int i=0; i<(G4int)radii.size(); i++){
		G4int cellID = i + 1;
		G4double count = 0.;
		G4double weighted = 0.;
		auto c = counts.find(cellID);
		auto w = weightedCounts.find(cellID);
		if(c != counts.end()) count = c->second;
		if(w != weightedCounts.end()) weighted = w->second;

		G4double currentImportance = (i + 1 < (G4int)importances.size()) ? importances[i + 1] : 1.;
		G4double suggestedImportance = (i + 1 < (G4int)recommended.size()) ? recommended[i + 1] : currentImportance;

		out << std::setw(5) << cellID
		       << " | " << std::setw(10) << radii[i]/cm
		       << " | " << std::setw(10) << count
		       << " | " << std::setw(13) << weighted
		       << " | " << std::setw(10) << currentImportance
		       << " | " << std::setw(12) << suggestedImportance
		       << G4endl;
	}

	out << " Recommended macro:" << G4endl
	       << "/tetcal/eyeImportanceRadii ";
	for(auto radius : radii) out << radius/cm << " ";
	out << G4endl << "/tetcal/eyeImportanceValues ";
	for(auto importance : recommended) out << importance << " ";
	out << G4endl
	       << "=======================================================================" << G4endl
	       << G4endl;
}

void EyeImportanceTuning::WriteTunedMacroAndData(const std::map<G4int, G4double>& counts,
                                                 const std::map<G4int, G4double>& weightedCounts,
                                                 G4int nEvents)
{
	if(!enabled || !tuneEnabled || radii.empty() || tuneInputMacro.empty() || tuneOutputMacro.empty()) return;

	auto recommended = RecommendImportances(weightedCounts, nEvents);
	G4double maxRelativeDifference = 0.;
	for(size_t i=0; i<recommended.size(); i++){
		G4double current = i < importances.size() ? importances[i] : 1.;
		maxRelativeDifference = std::max(maxRelativeDifference, RelativeDifference(current, recommended[i]));
	}
	tuneConverged = maxRelativeDifference <= tuneTolerance;

	std::ifstream ifs(tuneInputMacro);
	if(ifs.is_open()){
		G4String target = tuneOutputMacro;
		G4String tmpTarget = target + ".tmp";
		std::ofstream ofs(tmpTarget);
		G4bool replaced = false;
		G4bool inserted = false;
		G4String line;
		while(std::getline(ifs, line)){
			G4String command = FirstCommand(line);
			if(command == "/tetcal/eyeImportanceValues"){
				if(!replaced){
					ofs << ImportanceLine(recommended) << G4endl;
					replaced = true;
				}
				continue;
			}
			ofs << line << G4endl;
			if(!replaced && !inserted &&
			   (command == "/tetcal/eyeImportanceRadii" || command == "/tetcal/eyeImportance")){
				ofs << ImportanceLine(recommended) << G4endl;
				inserted = true;
			}
		}
		if(!replaced && !inserted){
			ofs << ImportanceLine(recommended) << G4endl;
		}
		ofs.close();
		if(ifs.good() || ifs.eof()){
			std::remove(target.c_str());
			std::rename(tmpTarget.c_str(), target.c_str());
		}
	}

	if(!tuneDataFile.empty()){
		std::ofstream data;
		if(tuneIteration <= 1){
			data.open(tuneDataFile);
			data << "iteration\tcell\tradius_cm\ttrack_enter\tweighted_enter\tcurrent_importance\tsuggested_importance\trelative_difference\tmax_relative_difference\tconverged" << G4endl;
		}
		else{
			data.open(tuneDataFile, std::ios::app);
		}

		for(G4int i=0; i<(G4int)radii.size(); i++){
			G4int cellID = i + 1;
			G4double count = 0.;
			G4double weighted = 0.;
			auto c = counts.find(cellID);
			auto w = weightedCounts.find(cellID);
			if(c != counts.end()) count = c->second;
			if(w != weightedCounts.end()) weighted = w->second;
			G4double current = (i + 1 < (G4int)importances.size()) ? importances[i + 1] : 1.;
			G4double suggested = (i + 1 < (G4int)recommended.size()) ? recommended[i + 1] : current;

			data << tuneIteration << "\t"
			     << cellID << "\t"
			     << std::setprecision(8) << radii[i]/cm << "\t"
			     << count << "\t"
			     << weighted << "\t"
			     << current << "\t"
			     << suggested << "\t"
			     << RelativeDifference(current, suggested) << "\t"
			     << maxRelativeDifference << "\t"
			     << (tuneConverged ? "yes" : "no") << G4endl;
		}
	}

	G4cout << "Eye importance tune iteration " << tuneIteration
	       << ": wrote " << tuneOutputMacro
	       << ", max relative change = " << maxRelativeDifference
	       << ", converged = " << (tuneConverged ? "yes" : "no")
	       << G4endl;
}
