#include "EyeImportanceParallelWorld.hh"

#include "TETModelImport.hh"

#include "G4GeometryCell.hh"
#include "G4IStore.hh"
#include "G4LogicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PVPlacement.hh"
#include "G4PSDirectionFlag.hh"
#include "G4PSTrackCounter.hh"
#include "G4SDManager.hh"
#include "G4SDParticleFilter.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

#include <algorithm>
#include <set>
#include <sstream>

EyeImportanceParallelWorld::EyeImportanceParallelWorld(const G4String& worldName,
                                                       TETModelImport* _tetData,
                                                       const std::vector<G4double>& radii,
                                                       const std::vector<G4double>& importances)
:G4VUserParallelWorld(worldName),
 tetData(_tetData),
 shellRadii(radii),
 shellImportances(importances),
 ghostWorld(nullptr)
{}

G4String EyeImportanceParallelWorld::GetCellName(G4int index) const
{
	std::ostringstream os;
	os << "eyeImportanceCell_" << index;
	return os.str();
}

G4bool EyeImportanceParallelWorld::ComputeEyeCenter(G4ThreeVector& center) const
{
	const std::set<G4int> leftLens = {6600, 6601};
	const std::set<G4int> rightLens = {6800, 6801};

	G4ThreeVector leftCenter;
	G4ThreeVector rightCenter;
	G4bool hasLeft = tetData->GetMaterialSetCenterOfMass(leftLens, leftCenter);
	G4bool hasRight = tetData->GetMaterialSetCenterOfMass(rightLens, rightCenter);

	if(hasLeft && hasRight){
		center = (leftCenter + rightCenter) * 0.5;
		return true;
	}
	if(hasLeft){
		center = leftCenter;
		return true;
	}
	if(hasRight){
		center = rightCenter;
		return true;
	}
	return false;
}

void EyeImportanceParallelWorld::Construct()
{
	ghostWorld = GetWorld();
	G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();
	logicalCells.push_back(worldLogical);

	if(!ComputeEyeCenter(eyeCenter)){
		G4Exception("EyeImportanceParallelWorld::Construct", "", JustWarning,
		            "Eye lens materials 6600/6601/6800/6801 were not found; eye importance geometry was not created.");
		return;
	}

	std::sort(shellRadii.begin(), shellRadii.end(), std::greater<G4double>());
	if(shellImportances.size() != shellRadii.size() + 1){
		G4Exception("EyeImportanceParallelWorld::Construct", "", FatalErrorInArgument,
		            "Eye importance needs one more importance value than the number of radii.");
	}

	G4VisAttributes* invisible = new G4VisAttributes();
	invisible->SetVisibility(false);

	for(G4int i=0; i<(G4int)shellRadii.size(); i++){
		G4double innerRadius = (i == (G4int)shellRadii.size() - 1) ? 0. : shellRadii[i + 1];
		G4double outerRadius = shellRadii[i];

		G4Sphere* sphere = new G4Sphere(GetCellName(i) + "_solid",
		                                innerRadius,
		                                outerRadius,
		                                0.*deg,
		                                360.*deg,
		                                0.*deg,
		                                180.*deg);
		G4LogicalVolume* logical = new G4LogicalVolume(sphere, nullptr, GetCellName(i) + "_logical");
		logical->SetVisAttributes(invisible);
		logicalCells.push_back(logical);

		G4VPhysicalVolume* physical = new G4PVPlacement(nullptr,
		                                                eyeCenter,
		                                                logical,
		                                                GetCellName(i),
		                                                worldLogical,
		                                                false,
		                                                i + 1);
		physicalCells.push_back(physical);
	}

	G4cout << "Eye importance center: "
	       << eyeCenter.x()/cm << " "
	       << eyeCenter.y()/cm << " "
	       << eyeCenter.z()/cm << " cm" << G4endl;
}

G4VIStore* EyeImportanceParallelWorld::CreateImportanceStore()
{
	G4IStore* istore = G4IStore::GetInstance(GetName());
	istore->AddImportanceGeometryCell(shellImportances.front(), *ghostWorld, 0);

	for(G4int i=0; i<(G4int)physicalCells.size(); i++){
		G4double importance = shellImportances[i + 1];
		istore->AddImportanceGeometryCell(importance, *physicalCells[i], i + 1);
		G4cout << "Eye importance cell " << i + 1
		       << " radius <= " << shellRadii[i]/cm << " cm"
		       << " importance " << importance << G4endl;
	}

	return istore;
}

void EyeImportanceParallelWorld::ConstructSD()
{
	if(logicalCells.empty()) return;

	G4SDManager* sdManager = G4SDManager::GetSDMpointer();
	G4MultiFunctionalDetector* detector = new G4MultiFunctionalDetector("EyeImportanceSD");
	sdManager->AddNewDetector(detector);

	G4SDParticleFilter* gammaFilter = new G4SDParticleFilter("EyeImportanceGammaFilter", "gamma");
	detector->SetFilter(gammaFilter);

	G4PSTrackCounter* trackEnter = new G4PSTrackCounter("TrackEnter", fCurrent_In);
	detector->RegisterPrimitive(trackEnter);

	G4PSTrackCounter* trackEnterWeighted = new G4PSTrackCounter("TrackEnterW", fCurrent_In);
	trackEnterWeighted->Weighted(true);
	detector->RegisterPrimitive(trackEnterWeighted);

	for(size_t i=1; i<logicalCells.size(); i++){
		SetSensitiveDetector(logicalCells[i], detector);
	}
}
