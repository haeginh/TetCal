/*
 * ExternalBeam.cc
 *
 *  Created on: Oct 18, 2019
 *      Author: hhg
 */

#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4RandomDirection.hh"
#include <set>
#include <algorithm>
#include <fstream>
#include "G4Box.hh"
#include "G4LogicalVolumeStore.hh"

#include "VOXModelImport.hh"
#include "SourceGenerator.hh"

InternalSource::InternalSource(VOXModelImport* _voxData)
:voxData(_voxData), isRBM(false)
{
    G4Box* phantomBox = (G4Box*) G4LogicalVolumeStore::GetInstance()->GetVolume("phantomLogical")->GetSolid();
    base = G4ThreeVector(-phantomBox->GetXHalfLength(),
                         -phantomBox->GetYHalfLength(),
                         -phantomBox->GetZHalfLength());
}

InternalSource::~InternalSource()
{}

void InternalSource::SetSource(std::vector<G4int> sources)
{
    std::set<G4int>    sourceSet(sources.begin(), sources.end());
    voxPick.clear(); spPick.clear(); isRBM = false;
    //Cout
    std::stringstream ss;
    ss<<"Set source organs for ";
    for(auto source:sourceSet) ss<<source<<" ";
    if(sources[0] == 0) {ss<<"(RBM)"; sourceSet.erase(0); isRBM = true;}

    //Extract source tet IDs
    if(!isRBM){
        for(G4int i=0;i<voxData->GetVoxelResolution(0);i++){
            for(G4int j=0;j<voxData->GetVoxelResolution(1);j++){
                for(G4int k=0;k<voxData->GetVoxelResolution(2);k++){
                    if(sourceSet.find(voxData->GetVoxelData(i, j, k))!= sourceSet.end())
                        voxPick[0].push_back(VOX(i, j, k));
                }
            }
        }
        ss<<" -> "<<voxPick[0].size()<<G4endl;
    }else{
        auto rbmRatio = voxData->GetRBMratio();
        for(auto sset:sourceSet){
            if(rbmRatio.find(sset)==rbmRatio.end()){
                G4cerr<<sset<<" is not included in RBMnBS file!"<<G4endl; exit(0);
            }
            spPick.push_back({rbmRatio[sset],sset});
        }
        std::sort(spPick.begin(), spPick.end());
        std::reverse(spPick.begin(), spPick.end());
        G4double previousW(0.);
        for(auto &sp:spPick) {
            sp.first += previousW;
            previousW = sp.first;
        }
        for(auto &sp:spPick) sp.first /= previousW;

        for(G4int i=0;i<voxData->GetVoxelResolution(0);i++){
            for(G4int j=0;j<voxData->GetVoxelResolution(1);j++){
                for(G4int k=0;k<voxData->GetVoxelResolution(2);k++){
                    if(sourceSet.find(voxData->GetVoxelData(i, j, k))!= sourceSet.end())
                        voxPick[voxData->GetVoxelData(i, j, k)].push_back(VOX(i, j, k));
                }
            }
        }
        G4int count(0);
        for(auto vp:voxPick) count += vp.second.size();
        ss<<" -> "<<count<<" ("<<voxPick.size()<<" regions)"<<G4endl;
    }
    G4cout<<ss.str();
}

void InternalSource::GetAprimaryPos(G4ThreeVector &position)
{
    G4int regionID(0);
    G4double rand = G4UniformRand();
    for(auto sp:spPick){
        if(rand>sp.first) continue;
        regionID = sp.second; break;
    }

    G4int randInt = floor(G4UniformRand() * (G4double)voxPick[regionID].size());
    position = RandomSamplingInAVoxel(voxPick[regionID][randInt]);
}

G4ThreeVector InternalSource::RandomSamplingInAVoxel(VOX vox){
    G4int i, j, k;
    std::tie(i, j, k) = vox;
    G4ThreeVector sampledPosition = G4ThreeVector((G4UniformRand()+i)*voxData->GetVoxelSize(0),
                                                  (G4UniformRand()+j)*voxData->GetVoxelSize(1),
                                                  (G4UniformRand()+k)*voxData->GetVoxelSize(2));
    return (sampledPosition + base);
}
