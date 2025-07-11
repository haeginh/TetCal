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


PrimaryGeneratorAction::PrimaryGeneratorAction(TETModelImport* _tetData)
:tetData(_tetData), fSourceGenerator(0), spectrumSource(false)
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

#include <fstream>
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