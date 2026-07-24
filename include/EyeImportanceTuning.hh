#ifndef EyeImportanceTuning_h
#define EyeImportanceTuning_h 1

#include "G4String.hh"
#include "globals.hh"

#include <map>
#include <ostream>
#include <vector>

class EyeImportanceTuning
{
public:
	static void Configure(G4bool enabled,
	                      const std::vector<G4double>& radii,
	                      const std::vector<G4double>& importances,
	                      const G4String& particleName);
	static void ConfigureTuneSession(G4bool enabled,
	                                 const G4String& inputMacro,
	                                 const G4String& tunedMacro,
	                                 const G4String& tuneData,
	                                 G4int iteration,
	                                 G4double tolerance);

	static G4bool IsEnabled();
	static G4bool IsTuneEnabled();
	static G4bool IsTuneConverged();
	static const std::vector<G4double>& GetRadii();
	static const std::vector<G4double>& GetImportances();
	static const G4String& GetParticleName();
	static const G4String& GetTunedMacroName();

	static std::vector<G4double> RecommendImportances(const std::map<G4int, G4double>& weightedCounts,
	                                                  G4int nEvents,
	                                                  G4double maxRatio = 4.);
	static void PrintSummary(std::ostream& out,
	                         const std::map<G4int, G4double>& counts,
	                         const std::map<G4int, G4double>& weightedCounts,
	                         G4int nEvents);
	static void WriteTunedMacroAndData(const std::map<G4int, G4double>& counts,
	                                   const std::map<G4int, G4double>& weightedCounts,
	                                   G4int nEvents);

private:
	static G4bool enabled;
	static G4bool tuneEnabled;
	static G4bool tuneConverged;
	static std::vector<G4double> radii;
	static std::vector<G4double> importances;
	static G4String particleName;
	static G4String tuneInputMacro;
	static G4String tuneOutputMacro;
	static G4String tuneDataFile;
	static G4int tuneIteration;
	static G4double tuneTolerance;
};

#endif
