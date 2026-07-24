#ifndef ShieldImportanceTuning_h
#define ShieldImportanceTuning_h 1

#include "G4String.hh"
#include "globals.hh"

#include <ostream>

class ShieldImportanceTuning
{
public:
	static void Configure(G4bool enabled,
	                      G4double importance,
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
	static G4double GetImportance();
	static const G4String& GetParticleName();
	static void ResetCounts();
	static void AddEnter(G4double weight);
	static void AddExit(G4double weight);
	static G4double GetEnter();
	static G4double GetWeightedEnter();
	static G4double GetExit();
	static G4double GetWeightedExit();

	static G4double RecommendImportance(G4double weightedEnter, G4double weightedExit);
	static void PrintSummary(std::ostream& out,
	                         G4double enter,
	                         G4double weightedEnter,
	                         G4double exit,
	                         G4double weightedExit);
	static void WriteTunedMacroAndData(G4double enter,
	                                   G4double weightedEnter,
	                                   G4double exit,
	                                   G4double weightedExit);

private:
	static G4bool enabled;
	static G4bool tuneEnabled;
	static G4bool tuneConverged;
	static G4double importance;
	static G4String particleName;
	static G4String tuneInputMacro;
	static G4String tuneOutputMacro;
	static G4String tuneDataFile;
	static G4int tuneIteration;
	static G4double tuneTolerance;
	static G4double enter;
	static G4double weightedEnter;
	static G4double exit;
	static G4double weightedExit;
};

#endif
