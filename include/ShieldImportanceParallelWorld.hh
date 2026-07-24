#ifndef ShieldImportanceParallelWorld_h
#define ShieldImportanceParallelWorld_h 1

#include "G4String.hh"
#include "G4VUserParallelWorld.hh"
#include "globals.hh"

#include <set>
#include <vector>

class G4LogicalVolume;
class G4VIStore;
class G4VPhysicalVolume;
class TETModelImport;

class ShieldImportanceParallelWorld : public G4VUserParallelWorld
{
public:
	ShieldImportanceParallelWorld(const G4String& worldName,
	                              TETModelImport* tetData,
	                              const std::set<G4int>& materialIDs,
	                              G4double importance);
	virtual ~ShieldImportanceParallelWorld() {}

	virtual void Construct();
	virtual void ConstructSD();

	G4VIStore* CreateImportanceStore();

private:
	TETModelImport* tetData;
	std::set<G4int> shieldMaterialIDs;
	G4double shieldImportance;
	std::vector<G4int> shieldCopyNumbers;
	G4VPhysicalVolume* ghostWorld;
	std::vector<G4VPhysicalVolume*> shieldPhysicals;
	std::vector<G4LogicalVolume*> shieldLogicals;
};

#endif
