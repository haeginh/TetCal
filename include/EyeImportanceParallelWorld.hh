#ifndef EyeImportanceParallelWorld_h
#define EyeImportanceParallelWorld_h 1

#include "G4GeometryCell.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4VUserParallelWorld.hh"
#include "globals.hh"

#include <vector>

class G4LogicalVolume;
class G4VIStore;
class G4VPhysicalVolume;
class TETModelImport;

class EyeImportanceParallelWorld : public G4VUserParallelWorld
{
public:
	EyeImportanceParallelWorld(const G4String& worldName,
	                           TETModelImport* tetData,
	                           const std::vector<G4double>& radii,
	                           const std::vector<G4double>& importances);
	virtual ~EyeImportanceParallelWorld() {}

	virtual void Construct();
	virtual void ConstructSD();

	G4VPhysicalVolume* GetWorldVolume() { return ghostWorld; }
	G4VIStore* CreateImportanceStore();

private:
	G4String GetCellName(G4int index) const;
	G4bool ComputeEyeCenter(G4ThreeVector& center) const;

	TETModelImport* tetData;
	std::vector<G4double> shellRadii;
	std::vector<G4double> shellImportances;
	std::vector<G4VPhysicalVolume*> physicalCells;
	std::vector<G4LogicalVolume*> logicalCells;
	G4VPhysicalVolume* ghostWorld;
	G4ThreeVector eyeCenter;
};

#endif
