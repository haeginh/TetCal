//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// TETPrimaryGeneratorAction.cc
// \file   MRCP_GEANT4/External/src/TETPrimaryGeneratorAction.cc
// \author Haegin Han
// \update
// \


#include "PrimaryGeneratorAction.hh"
#include "Randomize.hh"
#include "G4ParticleTable.hh"
#include "G4RandomDirection.hh"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>


PrimaryGeneratorAction::PrimaryGeneratorAction(TETModelImport* _tetData)
:tetData(_tetData), fSourceGenerator(0), fExternal(0), fInternal(0), fSurface(0),
 sourceName(""), spectrumSource(false), iaeaSpectrumSource(false)
{
	fParticleGun = new G4ParticleGun(1);
	fMessenger   = new PrimaryMessenger(this);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete fParticleGun;
	delete fMessenger;
    if(fExternal) delete fExternal;
    if(fInternal) delete fInternal;
    if(fSurface)  delete fSurface;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	if (iaeaSpectrumSource){
		GenerateIAEADecay(anEvent);
		return;
	}

	G4ThreeVector direction, position;
	fSourceGenerator->GetAprimaryPosDir(position, direction);
	fParticleGun->SetParticlePosition(position);
    fParticleGun->SetParticleMomentumDirection(direction);
	if (spectrumSource){
		fParticleGun->SetParticleEnergy(samplingE.lower_bound(G4UniformRand())->second);
		// G4cout<<samplingE.lower_bound(G4UniformRand())->second<<G4endl;
	}
	
	fParticleGun->GeneratePrimaryVertex(anEvent);
}

void PrimaryGeneratorAction::SetSpectrumSource(G4String inputFile)
{
	std::ifstream file(inputFile);
	if(!file.is_open()) {
		G4Exception("PrimaryGeneratorAction::SetSpectrumSource",
		            "FileNotFound", FatalException,
		            ("Cannot open file: " + inputFile).c_str());
		return;
	}

	spectrumSource = true;
	iaeaSpectrumSource = false;

	std::string::size_type idx = inputFile.rfind('.');
	if (idx != std::string::npos) {
		std::string extension = inputFile.substr(idx + 1);
		G4double pdfTot(0);
		std::vector<std::pair<G4double, G4double>> pdfVec;
		if (extension == "RAD" || extension == "rad") {
			if(RADcodes.empty()) {
				G4Exception("PrimaryGeneratorAction::SetSpectrumSource",
				            "NoRadCodes", FatalException,
				            "No RAD codes specified. Use /spec/RADcodes command.");
				return;
			}
			// Skip lines until "START RADIATION RECORDS" is found
			std::string line;
			while (std::getline(file, line)) {
				if (line.find("START RADIATION RECORDS") != std::string::npos) {
					break;
				}
			}

			int dummyInt;
			std::string type;
			G4double pdf, energy;
			while (file >> std::ws && !file.eof()) {
				file >> dummyInt >> pdf >> energy >> type;
				// 만약 type이 RADcodes에 포함되어 있지 않으면 진행
				if (std::find(RADcodes.begin(), RADcodes.end(), type) == RADcodes.end()) continue;
				// G4cout<<dummyInt<<" "<<pdf<<" "<<energy<<" "<<type<<G4endl;
				pdfVec.push_back(std::make_pair(energy*MeV, pdf));
				pdfTot += pdf;
			}
			samplingE.clear();
			G4double cdf(0);
			for (const auto& pair : pdfVec) {
				cdf += pair.second / pdfTot; // Cumulative distribution
				samplingE[cdf] = pair.first; // Store energy at CDF
			}
		}
		// 다른 파일은 다음에 작성 (아직 베타랑 오제전자를 자세히 해야할지 결정하지 못함)
		else{
			G4Exception("PrimaryGeneratorAction::SetSpectrumSource",
						"UnknownFileType", FatalException,
						("Unknown file extension: " + inputFile).c_str());
			return;
		}
	}
}

void PrimaryGeneratorAction::SetIAEASpectrumSource(G4String inputCommand)
{
	std::vector<G4String> args = SplitCommand(inputCommand);
	if(args.size()!=2){
		G4Exception("PrimaryGeneratorAction::SetIAEASpectrumSource",
					"WrongArgument", FatalException,
					"Use /spec/IAEA [directory] [nuclide], e.g. /spec/IAEA iaea_spectra I-131");
		return;
	}

	G4String directory = args[0];
	G4String nuclide = args[1];
	if(!directory.empty() && directory.back()!='/') directory += "/";

	iaeaDiscreteEmissions.clear();
	iaeaBetaIntervals.clear();
	iaeaBetaYield.clear();

	LoadIAEADiscrete(directory + nuclide + "_per_decay_discrete_aggregated.tsv");
	LoadIAEABeta(directory + nuclide + "_beta_continuum.tsv");

	if(iaeaDiscreteEmissions.empty() && iaeaBetaIntervals.empty()){
		G4Exception("PrimaryGeneratorAction::SetIAEASpectrumSource",
					"NoSpectrum", FatalException,
					("No usable IAEA spectrum records for " + nuclide).c_str());
		return;
	}

	spectrumSource = true;
	iaeaSpectrumSource = true;
	G4cout << "IAEA per-decay spectrum loaded for " << nuclide
		   << ": " << iaeaDiscreteEmissions.size() << " discrete records, "
		   << iaeaBetaIntervals.size() << " beta intervals" << G4endl;
}

std::vector<G4String> PrimaryGeneratorAction::SplitCommand(G4String input) const
{
	std::vector<G4String> tokens;
	G4String current;
	G4bool inQuote = false;
	char quoteChar = 0;
	for(size_t i=0; i<input.size(); i++){
		char c = input[i];
		if(inQuote){
			if(c==quoteChar) inQuote = false;
			else current += c;
		}
		else if(c=='"' || c=='\''){
			inQuote = true;
			quoteChar = c;
		}
		else if(std::isspace(static_cast<unsigned char>(c))){
			if(!current.empty()){
				tokens.push_back(current);
				current.clear();
			}
		}
		else current += c;
	}
	if(!current.empty()) tokens.push_back(current);
	return tokens;
}

G4ParticleDefinition* PrimaryGeneratorAction::GetParticleDefinition(G4String particleName) const
{
	if(particleName=="gamma" || particleName=="xray") particleName = "gamma";
	else if(particleName=="electron" || particleName=="beta-") particleName = "e-";
	else if(particleName=="beta+") particleName = "e+";
	else if(particleName=="alpha") particleName = "alpha";
	else return 0;

	return G4ParticleTable::GetParticleTable()->FindParticle(particleName);
}

void PrimaryGeneratorAction::LoadIAEADiscrete(G4String inputFile)
{
	std::ifstream file(inputFile);
	if(!file.is_open()){
		G4Exception("PrimaryGeneratorAction::LoadIAEADiscrete",
					"FileNotFound", FatalException,
					("Cannot open file: " + inputFile).c_str());
		return;
	}

	G4String line;
	if(!std::getline(file, line)) return;
	std::vector<G4String> headers = SplitCommand(line);
	std::map<G4String, size_t> column;
	for(size_t i=0; i<headers.size(); i++) column[headers[i]] = i;

	while(std::getline(file, line)){
		if(line.empty()) continue;
		std::vector<G4String> fields;
		std::stringstream ss(line);
		G4String field;
		while(std::getline(ss, field, '\t')) fields.push_back(field);
		if(fields.size()<=column["yield_per_decay"] || fields.size()<=column["energy_MeV"]) continue;

		G4ParticleDefinition* particle = GetParticleDefinition(fields[column["particle"]]);
		if(!particle) continue;

		G4double yield = std::atof(fields[column["yield_per_decay"]].c_str());
		G4double energy = std::atof(fields[column["energy_MeV"]].c_str())*MeV;
		if(yield<=0. || energy<0.) continue;

		iaeaDiscreteEmissions.push_back({particle, energy, yield});
	}
}

void PrimaryGeneratorAction::LoadIAEABeta(G4String inputFile)
{
	std::ifstream file(inputFile);
	if(!file.is_open()){
		G4Exception("PrimaryGeneratorAction::LoadIAEABeta",
					"FileNotFound", FatalException,
					("Cannot open file: " + inputFile).c_str());
		return;
	}

	G4String line;
	if(!std::getline(file, line)) return;
	std::vector<G4String> headers = SplitCommand(line);
	std::map<G4String, size_t> column;
	for(size_t i=0; i<headers.size(); i++) column[headers[i]] = i;

	std::map<G4ParticleDefinition*, std::vector<std::pair<G4double, G4double>>> points;
	while(std::getline(file, line)){
		if(line.empty()) continue;
		std::vector<G4String> fields;
		std::stringstream ss(line);
		G4String field;
		while(std::getline(ss, field, '\t')) fields.push_back(field);
		if(fields.size()<=column["particle"] || fields.size()<=column["bin_energy_MeV"] ||
		   fields.size()<=column["dn_dE_per_decay_per_MeV"]) continue;

		G4ParticleDefinition* particle = GetParticleDefinition(fields[column["particle"]]);
		if(!particle) continue;

		G4double energy = std::atof(fields[column["bin_energy_MeV"]].c_str())*MeV;
		G4double dnDe = std::atof(fields[column["dn_dE_per_decay_per_MeV"]].c_str())/MeV;
		points[particle].push_back({energy, dnDe});
	}

	for(auto& entry : points){
		auto& vec = entry.second;
		std::sort(vec.begin(), vec.end());
		for(size_t i=0; i+1<vec.size(); i++){
			G4double e0 = vec[i].first;
			G4double e1 = vec[i+1].first;
			G4double y0 = vec[i].second;
			G4double y1 = vec[i+1].second;
			G4double weight = 0.5*(y0+y1)*(e1-e0);
			if(weight<=0.) continue;
			iaeaBetaIntervals.push_back({entry.first, e0, e1, weight});
			iaeaBetaYield[entry.first] += weight;
		}
	}
}

void PrimaryGeneratorAction::GenerateIAEADecay(G4Event* anEvent)
{
	if(!fSourceGenerator){
		G4Exception("PrimaryGeneratorAction::GenerateIAEADecay",
					"NoSource", FatalException,
					"Set /internal/source, /internal/surface, or /external/dir before /run/beamOn.");
		return;
	}

	G4ThreeVector position, direction;
	fSourceGenerator->GetAprimaryPosDir(position, direction);

	for(const auto& emission : iaeaDiscreteEmissions){
		G4int multiplicity = SampleMultiplicity(emission.yield);
		for(G4int i=0; i<multiplicity; i++)
			GenerateOnePrimary(anEvent, position, emission.particle, emission.energy);
	}

	for(const auto& by : iaeaBetaYield){
		G4int multiplicity = SampleMultiplicity(by.second);
		for(G4int i=0; i<multiplicity; i++)
			GenerateOnePrimary(anEvent, position, by.first, SampleBetaEnergy(by.first));
	}
}

void PrimaryGeneratorAction::GenerateOnePrimary(G4Event* anEvent,
                                                const G4ThreeVector& position,
                                                G4ParticleDefinition* particle,
                                                G4double energy)
{
	fParticleGun->SetParticleDefinition(particle);
	fParticleGun->SetParticlePosition(position);
	fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
	fParticleGun->SetParticleEnergy(energy);
	fParticleGun->GeneratePrimaryVertex(anEvent);
}

G4int PrimaryGeneratorAction::SampleMultiplicity(G4double yield) const
{
	if(yield<=0.) return 0;
	G4int multiplicity = G4int(std::floor(yield));
	G4double residual = yield - multiplicity;
	if(G4UniformRand()<residual) multiplicity++;
	return multiplicity;
}

G4double PrimaryGeneratorAction::SampleBetaEnergy(G4ParticleDefinition* particle) const
{
	auto itYield = iaeaBetaYield.find(particle);
	if(itYield==iaeaBetaYield.end() || itYield->second<=0.) return 0.;

	G4double rand = G4UniformRand()*itYield->second;
	G4double cdf = 0.;
	for(const auto& interval : iaeaBetaIntervals){
		if(interval.particle!=particle) continue;
		cdf += interval.weight;
		if(rand<=cdf)
			return interval.energy0 + G4UniformRand()*(interval.energy1-interval.energy0);
	}
	return 0.;
}
