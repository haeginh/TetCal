#ifndef INCLUDE_IMPORTVOXELPHANTOM_HH_
#define INCLUDE_IMPORTVOXELPHANTOM_HH_

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>

#include "G4ios.hh"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Element.hh"

class VOXModelImport {
public:
    VOXModelImport(G4String _phantomFile);
    virtual ~VOXModelImport();

    void ImportPhantomInfo(G4String);
    void ImportPhantomVoxelData(G4String);
    void ImportPhantomVoxelMaterial(G4String);
    void ImportPhantomVoxelVolume();

    G4String GetInputFilename() {return Filename;}
    G4int GetVoxelResolution(G4int idx) {return voxelResolution[idx];}
    G4double GetVoxelSize(G4int idx) {return voxelSize[idx];}
    G4Material* GetVoxelMaterial(G4int idx) {return materialMap[idx];}
    G4double GetOrganVolume(G4int idx) {return organVolume[idx];}
    G4double GetOrganMass(G4int idx) {return organVolume[idx] * materialMap[idx]->GetDensity();}
    G4int GetMaterialIndex(G4int idx) {return organID[idx];}
    std::map<G4int, G4double> GetMassMap() {return massMap;}
    G4int GetVoxelMaterialSize() {return (G4int) materialMap.size();}
    G4int GetVoxelData(G4int idx, G4int idy, G4int idz) {return voxelData[idx][idy][idz];}
    G4String GetOrganName(G4int idx) {return organNameMap[idx];}
    std::map<G4int, G4double> GetRBMratio()  { return rbmRatio;}

private:
    std::map<G4int, G4double> organVolume;
    std::map<G4int, G4int> numVoxel;
    std::map<G4int, std::vector<std::pair<G4int, G4double> > > materialIndexMap;
    std::vector<G4int> organID;
    std::map<G4int, G4Material* > materialMap;
    std::map<G4int, G4double> densityMap;
    std::map<G4int, G4double> massMap;
    std::map<G4int, G4String> organNameMap;

    G4ThreeVector voxelSize;
    std::vector<G4int> voxelResolution;
    G4ThreeVector phantomSize;

    G4int*** voxelData;

    G4String Filename;

    void RBMBSRead(G4String);

    std::map<G4int, G4double>  rbmRatio;
};

#endif /* INCLUDE_IMPORTVOXELPHANTOM_HH_ */
