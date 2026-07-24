#include "ShieldImportanceParallelWorld.hh"

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
#include "G4SystemOfUnits.hh"
#include "G4Tet.hh"
#include "G4VisAttributes.hh"

ShieldImportanceParallelWorld::ShieldImportanceParallelWorld(const G4String& worldName,
                                                             TETModelImport* _tetData,
                                                             const std::set<G4int>& materialIDs,
                                                             G4double importance)
: G4VUserParallelWorld(worldName),
  tetData(_tetData),
  shieldMaterialIDs(materialIDs),
  shieldImportance(importance),
  ghostWorld(nullptr)
{}

void ShieldImportanceParallelWorld::Construct()
{
	ghostWorld = GetWorld();
	G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();

	shieldCopyNumbers = tetData->GetCopyNumbersForMaterials(shieldMaterialIDs);
	if(shieldCopyNumbers.empty()){
		G4Exception("ShieldImportanceParallelWorld::Construct", "", JustWarning,
		            "No lead/lead_glass tetrahedra were found; shield importance geometry was not created.");
		return;
	}

	G4VisAttributes* invisible = new G4VisAttributes();
	invisible->SetVisibility(false);

	for(G4int i=0; i<(G4int)shieldCopyNumbers.size(); i++){
		G4LogicalVolume* logical = new G4LogicalVolume(tetData->GetTetrahedron(shieldCopyNumbers[i]),
		                                               nullptr,
		                                               "ShieldImportanceLogical");
		logical->SetVisAttributes(invisible);
		shieldLogicals.push_back(logical);

		G4VPhysicalVolume* physical = new G4PVPlacement(nullptr,
		                                                G4ThreeVector(),
		                                                logical,
		                                                "ShieldImportancePhysical",
		                                                worldLogical,
		                                                false,
		                                                i);
		shieldPhysicals.push_back(physical);
	}

	G4cout << "Shield importance cells: " << shieldCopyNumbers.size()
	       << " tetrahedra for materials";
	for(auto id : shieldMaterialIDs) G4cout << " " << id;
	G4cout << G4endl;
}

G4VIStore* ShieldImportanceParallelWorld::CreateImportanceStore()
{
	G4IStore* istore = G4IStore::GetInstance(GetName());
	istore->AddImportanceGeometryCell(1., *ghostWorld, 0);

	for(G4int i=0; i<(G4int)shieldPhysicals.size(); i++){
		istore->AddImportanceGeometryCell(shieldImportance, *shieldPhysicals[i], 0);
	}
	G4cout << "Shield importance value: " << shieldImportance << G4endl;

	return istore;
}

void ShieldImportanceParallelWorld::ConstructSD()
{
	if(shieldLogicals.empty()) return;

	G4SDManager* sdManager = G4SDManager::GetSDMpointer();
	G4MultiFunctionalDetector* detector = new G4MultiFunctionalDetector("ShieldImportanceSD");
	sdManager->AddNewDetector(detector);

	G4SDParticleFilter* gammaFilter = new G4SDParticleFilter("ShieldImportanceGammaFilter", "gamma");
	detector->SetFilter(gammaFilter);

	G4PSTrackCounter* trackEnter = new G4PSTrackCounter("TrackEnter", fCurrent_In);
	detector->RegisterPrimitive(trackEnter);

	G4PSTrackCounter* trackEnterWeighted = new G4PSTrackCounter("TrackEnterW", fCurrent_In);
	trackEnterWeighted->Weighted(true);
	detector->RegisterPrimitive(trackEnterWeighted);

	G4PSTrackCounter* trackExit = new G4PSTrackCounter("TrackExit", fCurrent_Out);
	detector->RegisterPrimitive(trackExit);

	G4PSTrackCounter* trackExitWeighted = new G4PSTrackCounter("TrackExitW", fCurrent_Out);
	trackExitWeighted->Weighted(true);
	detector->RegisterPrimitive(trackExitWeighted);

	for(auto logical : shieldLogicals) SetSensitiveDetector(logical, detector);
}
