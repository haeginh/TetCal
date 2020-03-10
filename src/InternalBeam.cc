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

#include "../include/VOXModelImport.hh"
#include "SourceGenerator.hh"

InternalSource::InternalSource(VOXModelImport* _voxData)
:voxData(_voxData)
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

	//Cout
	std::stringstream ss;
	ss<<"Set source organs for ";
	for(auto source:sourceSet) ss<<source<<" ";

	//Extract source tet IDs
	for(G4int i=0;i<voxData->GetVoxelResolution(0);i++){
		for(G4int j=0;j<voxData->GetVoxelResolution(1);j++){
			for(G4int k=0;k<voxData->GetVoxelResolution(2);k++){
				if(sourceSet.find(voxData->GetVoxelData(i, j, k))!= sourceSet.end())
					voxPick.push_back(VOX(i, j, k));
			}
		}
	}
	ss<<" -> "<<voxPick.size()<<G4endl;
	G4cout<<ss.str();
}

void InternalSource::GetAprimaryPosDir(G4ThreeVector &position, G4ThreeVector &direction)
{
	G4int randInx = floor(G4UniformRand() * (G4double)voxPick.size());
	position = RandomSamplingInAVoxel(voxPick[randInx]);
	direction = G4RandomDirection();
}

G4ThreeVector InternalSource::RandomSamplingInAVoxel(VOX vox){
	G4int i, j, k;
	std::tie(i, j, k) = vox;
	G4ThreeVector sampledPosition = G4ThreeVector((G4UniformRand()+i)*voxData->GetVoxelSize(0),
												  (G4UniformRand()+j)*voxData->GetVoxelSize(1),
												  (G4UniformRand()+k)*voxData->GetVoxelSize(2));
	return (sampledPosition + base);
}


