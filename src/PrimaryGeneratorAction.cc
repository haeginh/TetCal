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


#include "PrimaryGeneratorAction.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
#include <fstream>

PrimaryGeneratorAction::PrimaryGeneratorAction(G4double _centerZ)
:centerZ(_centerZ)
{
	fParticleGun = new G4ParticleGun(1);
	fMessenger   = new PrimaryMessenger(this);
	fParticleGun->SetParticleDefinition(G4Gamma::GammaDefinition());
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete fParticleGun;
	delete fMessenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	G4double rand = G4UniformRand();
	G4double rand_energy = cdf.rbegin()->second;
	if (rand < 1)
	{
		for (auto itr : cdf)
		{
			if (rand > itr.first) continue;
			rand_energy = itr.second;
			break;
		}
	}
	fParticleGun->SetParticleEnergy(rand_energy);

	G4ThreeVector direction, position;
	position.setZ(lowerBound+(upperBound-lowerBound)*G4UniformRand());
	G4double theta = G4UniformRand()*2*M_PI;
	position.setX(radius*cos(theta));
	position.setY(radius*sin(theta));
	fParticleGun->SetParticlePosition(position);

	G4double theta1 = -angle*0.5 + angle*G4UniformRand();
	direction.setX(-cos(theta+theta1));
	direction.setY(-sin(theta+theta1));
	fParticleGun->SetParticleMomentumDirection(direction);

	fParticleGun->GeneratePrimaryVertex(anEvent);
}

void PrimaryGeneratorAction::SetPeakEnergy(G4int _kVp)
{
	using namespace std;
	kVp = _kVp;
	G4String fileName(to_string(kVp) + ".spec");
	G4String spectra(specDir + "/CT" + fileName);

	G4cout << "Read x-ray spectra: " << spectra << G4endl;
	ifstream ifs(spectra);

	vector<pair<G4double, G4double>> pdf; //there could be duplicated probabilities, which may cause error in 'map'
	if (!ifs.is_open())
	{
		G4cerr << "X-ray spectra file was not opened" << G4endl;
		exit(1);
	}

	specSum = 0;
	G4String dump;
	while (getline(ifs, dump))
	{
		stringstream ss(dump);
		ss >> dump;
		if (dump == "Energy[keV]")
		{
			while (getline(ifs, dump))
			{
				G4double energy, intensity;
				stringstream ss2(dump);
				ss2 >> energy >> intensity;
				specSum += intensity;
				pdf.push_back(make_pair(intensity, energy * keV));
			}
		}
	}
	ifs.close();

	sort(pdf.begin(), pdf.end(), greater<std::pair<G4double, G4double>>());

	cdf.clear();
	G4double sumProb(0);
	for (auto itr : pdf)
	{
		sumProb += itr.first/specSum;
		cdf[sumProb] = itr.second;
	}
}
