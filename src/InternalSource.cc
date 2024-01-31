/*
 * ExternalBeam.cc
 *
 *  Created on: Oct 18, 2019
 *      Author: hhg
 */

#include "SeedParallel.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4RandomDirection.hh"
#include <set>
#include <algorithm>
#include <fstream>
#include "SourceGenerator.hh"

InternalSource::InternalSource(SeedParallel* _seedParallel)
:seedParallel(_seedParallel)
{}

InternalSource::~InternalSource()
{}

void InternalSource::SetSource(std::vector<G4int> sources)
{
	std::set<G4int>    sourceSet(sources.begin(), sources.end());
	tetPick.clear();
	
	//Cout
	std::stringstream ss;
	ss<<"Set source organs for "<<G4endl;
	for(auto source:sourceSet) ss<<source<<" ";

	//Extract source tet IDs
	for(G4int i=0;i<seedParallel->GetNumOfTet();i++){
		if(sourceSet.find(seedParallel->GetId(i)) != sourceSet.end())
			tetPick.push_back(VOLPICK(seedParallel->GetTet(i)->GetCubicVolume(), i));
	}
    ss<<" -> "<<tetPick.size()<<G4endl;
	if(tetPick.size()==0){
		G4Exception("InternalSource::SetSource","",FatalErrorInArgument,
				G4String("       Wrong source ID wad defined" ).c_str());
	}

	//Arrange volumes
	std::sort(tetPick.begin(), tetPick.end());
	std::reverse(tetPick.begin(), tetPick.end());

	G4double previousVol(0.);
	for(auto &tp:tetPick) {
		tp.first += previousVol;
		previousVol = tp.first;
	}

	for(auto &tp:tetPick) tp.first /= previousVol;
	G4cout<<ss.str();
}

void InternalSource::GetAprimaryPosDir(G4ThreeVector &position, G4ThreeVector &direction)
{
	G4double rand = G4UniformRand();
	for(auto tp:tetPick){
		if(rand>tp.first) continue;
		position = RandomSamplingInTet(seedParallel->GetTet(tp.second)); break;
	}
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
	return SampledPosition + seedParallel->GetCenter();
}


