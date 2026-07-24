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
#include <set>
#include <algorithm>
#include <fstream>
#include <cctype>
#include "SourceGenerator.hh"

InternalSource::InternalSource(TETModelImport* _tetData)
:tetData(_tetData)
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

	std::vector<G4int> copyNumbers;
	for(G4int i=0;i<tetData->GetNumTetrahedron();i++){
		if(sourceSet.find(tetData->GetMaterialIndex(i)) != sourceSet.end())
			copyNumbers.push_back(i);
	}
	BuildSourceFromCopyNumbers(copyNumbers, ss.str());
}

void InternalSource::SetSourceElements(std::vector<G4int> elementIDs, std::vector<G4int> materialIDs)
{
	std::set<G4int> elementSet(elementIDs.begin(), elementIDs.end());
	std::set<G4int> materialSet(materialIDs.begin(), materialIDs.end());
	tetPick.clear();

	std::stringstream ss;
	ss<<"Set source element IDs for "<<G4endl;
	ss<<elementSet.size()<<" element IDs";
	if(!materialSet.empty()){
		ss<<" with material filter ";
		for(auto materialID : materialSet) ss<<materialID<<" ";
	}

	BuildSourceFromCopyNumbers(tetData->GetCopyNumbersForElementIDs(elementSet, materialSet), ss.str());
}

void InternalSource::SetSourceElementFile(G4String sourceFile, std::vector<G4int> materialIDs)
{
	SetSourceElements(ReadElementIDFile(sourceFile), materialIDs);
}

void InternalSource::BuildSourceFromCopyNumbers(const std::vector<G4int>& copyNumbers, const G4String& sourceDescription)
{
	for(auto copyNo : copyNumbers){
		if(copyNo < 0 || copyNo >= tetData->GetNumTetrahedron()) continue;
		tetPick.push_back(VOLPICK(tetData->GetTetrahedron(copyNo)->GetCubicVolume(), copyNo));
	}

	std::stringstream ss;
	ss<<sourceDescription<<" -> "<<tetPick.size()<<G4endl;
	if(tetPick.size()==0){
		G4Exception("InternalSource::SetSource","",FatalErrorInArgument,
				G4String("       No source tetrahedron was selected" ).c_str());
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

std::vector<G4int> InternalSource::ReadElementIDFile(G4String sourceFile)
{
	std::ifstream ifs(sourceFile);
	if(!ifs.is_open()){
		G4Exception("InternalSource::ReadElementIDFile","",FatalErrorInArgument,
				G4String("       There is no source element file: " + sourceFile).c_str());
	}

	std::vector<G4int> elementIDs;
	G4String line;
	while(std::getline(ifs, line)){
		auto commentPos = line.find('#');
		if(commentPos != std::string::npos) line = line.substr(0, commentPos);

		line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](unsigned char ch){ return !std::isspace(ch); }));
		line.erase(std::find_if(line.rbegin(), line.rend(), [](unsigned char ch){ return !std::isspace(ch); }).base(), line.end());
		if(line.empty() || line.front() == '[') continue;

		auto dash = line.find('-');
		if(dash == std::string::npos){
			elementIDs.push_back(std::atoi(line.c_str()));
			continue;
		}

		G4int firstID = std::atoi(line.substr(0, dash).c_str());
		G4int lastID = std::atoi(line.substr(dash + 1).c_str());
		if(lastID < firstID) std::swap(firstID, lastID);
		for(G4int elementID = firstID; elementID <= lastID; elementID++) elementIDs.push_back(elementID);
	}
	return elementIDs;
}

void InternalSource::GetAprimaryPosDir(G4ThreeVector &position, G4ThreeVector &direction)
{
	G4double rand = G4UniformRand();
	for(auto tp:tetPick){
		if(rand>tp.first) continue;
		position = RandomSamplingInTet(tetData->GetTetrahedron(tp.second)); break;
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
	return SampledPosition;
}

