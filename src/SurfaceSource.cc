/*
 * SurfaceSource.cc
 *
 *  Created on: Mar 20, 2020
 *      Author: hhg
 */

#include "TETModelImport.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4RandomDirection.hh"
#include <set>
#include <algorithm>
#include <fstream>
#include <cmath>
#include "SourceGenerator.hh"

SurfaceSource::SurfaceSource(TETModelImport* _tetData)
:tetData(_tetData)
{}

SurfaceSource::~SurfaceSource()
{}

void SurfaceSource::SetSource(std::vector<G4int> sources)
{
    facePick.clear();
    for(G4int source:sources){
         std::vector<std::tuple<G4int, G4int, G4int>> facePool;
        auto elements = tetData->GetElements(source);
        for(auto ele:elements){
            facePool.push_back({ele[0], ele[1], ele[2]});
            facePool.push_back({ele[0], ele[1], ele[3]});
            facePool.push_back({ele[0], ele[2], ele[3]});
            facePool.push_back({ele[1], ele[2], ele[3]});
        }

        sort(facePool.begin(), facePool.end());
        for(size_t i=0;i<facePool.size()-1;i++){
            if(facePool[i]==facePool[i+1]){
                i++; continue;
            }
            facePick.push_back({CalculateTriangleArea(tetData->GetAVertex(std::get<0>(facePool[i])),
                                                      tetData->GetAVertex(std::get<1>(facePool[i])),
                                                      tetData->GetAVertex(std::get<2>(facePool[i]))), facePool[i]});
            if(i==facePool.size()-2){
                facePick.push_back({CalculateTriangleArea(tetData->GetAVertex(std::get<0>(facePool[i+1])),
                                                          tetData->GetAVertex(std::get<1>(facePool[i+1])),
                                                          tetData->GetAVertex(std::get<2>(facePool[i+1]))), facePool[i+1]});
            }
        }
    }
    if(facePick.empty()){
  		G4Exception("SurfaceSource::SetSource","",FatalErrorInArgument,
		G4String("       Wrong source ID wad defined" ).c_str());
    }

    std::sort(facePick.begin(), facePick.end());
    std::reverse(facePick.begin(), facePick.end());

    G4double previousVol(0.);
    for(auto &fp:facePick) {
        fp.first += previousVol;
        previousVol = fp.first;
    }

    for(auto &fp:facePick) fp.first /= previousVol;
    std::stringstream ss;
    ss<<"Surface Source: ";
    for(auto source:sources) ss<<source<<" ";
    ss<<"-> "<<facePick.size()<<" faces"<<G4endl;
    G4cout<<ss.str();
}

void SurfaceSource::GetAprimaryPosDir(G4ThreeVector &position, G4ThreeVector &direction)
{
    G4double rand = G4UniformRand();
    for(auto fp:facePick){
        if(rand>fp.first) continue;
        position = RandomSamplingInTriangle(fp.second); break;
    }
    direction = G4RandomDirection();
}

G4double SurfaceSource::CalculateTriangleArea(G4ThreeVector a, G4ThreeVector b, G4ThreeVector c){
    return 0.5*(b-a).cross(c-a).mag();
}

G4ThreeVector SurfaceSource::RandomSamplingInTriangle(std::tuple<G4int, G4int, G4int> triangle){

    G4double rand1 = G4UniformRand();
    G4double rand2 = G4UniformRand();

    if (rand1+rand2>1.0){
        rand1 = 1.0 - rand1;
        rand2 = 1.0 - rand2;
    }

    G4ThreeVector vecA = tetData->GetAVertex(std::get<1>(triangle)) - tetData->GetAVertex(std::get<0>(triangle));
    G4ThreeVector vecB = tetData->GetAVertex(std::get<2>(triangle)) - tetData->GetAVertex(std::get<0>(triangle));

    return tetData->GetAVertex(std::get<0>(triangle)) + rand1*vecA + rand2*vecB;
 }


