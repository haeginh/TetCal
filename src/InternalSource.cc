/*
 * ExternalBeam.cc
 *
 *  Created on: Oct 18, 2019
 *      Author: hhg
 */

#include "TETModelImport.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4RandomDirection.hh"
#include <algorithm>
#include <exception>
#include <unordered_set>
#include "SourceGenerator.hh"
#include "SourceDistribution.hh"

InternalSource::InternalSource(TETModelImport* _tetData)
:useFractions(false), tetData(_tetData)
{}

InternalSource::~InternalSource()
{}

void InternalSource::SetSource(std::vector<G4int> sources)
{
	sourceIDs = sources;
	sourceFractions.clear();
	useFractions = false;
	tetPick.clear();
	G4cout << "Configured internal source organs: ";
	for(auto source : sourceIDs) G4cout << source << " ";
	G4cout << G4endl;
}

void InternalSource::SetFractions(std::vector<G4double> fractions)
{
	if(sourceIDs.empty()){
		G4Exception("InternalSource::SetFractions", "", FatalErrorInArgument,
				"/internal/source must be defined before /internal/fraction");
		return;
	}
	sourceFractions = fractions;
	useFractions = true;
	BuildTetPick();
	G4cout << "Set fraction-weighted source organs with "
	       << sourceFractions.size() << " fractions -> "
	       << tetPick.size() << " tetrahedra" << G4endl;
}

void InternalSource::BuildTetPick()
{
	std::unordered_set<G4int> requested(sourceIDs.begin(), sourceIDs.end());
	std::vector<tetcal::SourceTet> tetrahedra;
	for(G4int i = 0; i < tetData->GetNumTetrahedron(); ++i){
		if(requested.find(tetData->GetMaterialIndex(i)) == requested.end()) continue;
		tetrahedra.push_back(tetcal::SourceTet{
			tetData->GetTetrahedron(i)->GetCubicVolume(),
			tetData->GetMaterialIndex(i),
			i});
	}

	try {
		std::vector<int> ids(sourceIDs.begin(), sourceIDs.end());
		if(!useFractions){
			tetPick = tetcal::BuildVolumeWeightedCDF(tetrahedra, ids);
		} else {
			std::vector<double> fractions(sourceFractions.begin(), sourceFractions.end());
			tetPick = tetcal::BuildFractionWeightedCDF(tetrahedra, ids, fractions);
		}
	} catch(const std::exception& error) {
		tetPick.clear();
		G4Exception("InternalSource::BuildTetPick", "", FatalErrorInArgument,
				error.what());
	}
	G4cout << (!useFractions ? "Built volume-weighted source over "
	                        : "Built fraction-weighted source over ")
	       << tetPick.size() << " tetrahedra" << G4endl;
}

void InternalSource::GetAprimaryPosDir(G4ThreeVector &position, G4ThreeVector &direction)
{
	if(tetPick.empty()) BuildTetPick();
	G4double rand = G4UniformRand();
	std::vector<VOLPICK>::const_iterator selected = std::lower_bound(
		tetPick.begin(), tetPick.end(), rand,
		[](const VOLPICK& entry, G4double value) {return entry.first < value;});
	if(selected == tetPick.end()) selected = tetPick.end() - 1;
	position = RandomSamplingInTet(tetData->GetTetrahedron(selected->second));
	direction = G4RandomDirection();
}

G4ThreeVector InternalSource::RandomSamplingInTet(G4Tet* tet){

	G4double varS = G4UniformRand();
	G4double varT = G4UniformRand();
	G4double varU = G4UniformRand();

	if (varS+varT>1.0){

		varS = 1.0 - varS;
		varT = 1.0 - varT;

	}
	if (varT+varU>1.0){

		double tmp = varU;
		varU = 1.0 - varS - varT;
		varT = 1.0 -tmp;
	} else if (varS+varT+varU>1.0){

		double tmp = varU;
		varU = varS + varT + varU - 1.0;
		varS = 1 - varT - tmp;
	}

	double a = 1 - varS - varT - varU;
	
	G4ThreeVector v1, v2, v3, v4;
	tet->GetVertices(v1, v2, v3, v4);
	G4ThreeVector SampledPosition = a*v1+varS*v2+varT*v3+varU*v4;
	return SampledPosition;
}
