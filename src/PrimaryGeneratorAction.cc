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
// The code was written by :
//	*Yeon Soo Yeom, yeonsoo.yeom@nih.gov
//	*Choonsik Lee, choonsik.lee@nih.gov
//
// Radiation Epidemiology Branch, DCEG/NCI/NHI
// 9609 Medical Center Dr, Rociville, MD 20850
// Tel: +1-240-276-532
// ********************************************************************

#include "PrimaryGeneratorAction.hh"

#include "G4GeneralParticleSource.hh"
#include "G4Event.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include <fstream>
#include "stdlib.h"
#include "time.h"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(G4String _inp)
:inp_name (_inp)
{
	GetSourceInfo(inp_name);
	particleGun = new G4ParticleGun();
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle;
	particle = particleTable->FindParticle("gamma");
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(PosX,PosY,PosZ));
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete particleGun;
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{	SourceSampling();
	particleGun->SetParticleMomentumDirection(direction);
	//particleGun->SetParticleEnergy(0.05*MeV);
	particleGun->GeneratePrimaryVertex(anEvent);
}

void PrimaryGeneratorAction::SourceSampling(){
	//Sampling Energy
	// G4double rand1= 1.0*rand()/RAND_MAX*CDF[CDF.size()-1];
	
	// if(CDF[0] >= rand1) {
	// energy= (1.0*rand()/RAND_MAX)*(PDF_Up[0].second-PDF_Dw[0].second) + PDF_Dw[0].second;	
	// }
	
	
	// for (int i=CDF.size()-1; i>0; i--){
	// 	if(CDF[i] >= rand1 && CDF[i-1] < rand1){
	// 	energy= (1.0*rand()/RAND_MAX)*(PDF_Up[i].second-PDF_Dw[i].second) + PDF_Dw[i].second;
	// 	break;		
	// 	}
	// }
	
	// G4double rand = G4UniformRand();
	// auto lb = CDF.lower_bound(rand); 
	// energy = (prev(lb))->second + (lb->second - prev(lb)->second)/(lb->first - prev(lb)->second)*rand;
	// //Sampling Direction
	
	G4double rand2= G4UniformRand();
	G4double theta = 2.0*pi*G4UniformRand();
	//G4double phi = acos(1-(1-ConeAngle)*(1.0*rand()/RAND_MAX));
	G4double phi = acos((ConeAngle-1.0)*rand2+1.0);
	//G4double phi = 0.0;
	G4ThreeVector V3(sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi));
	//G4ThreeVector V3(0, 0, 1);
	direction = V3.rotate(RotAngle, RotAxis);
	
	//cout << direction.angle(G4ThreeVector(VecX,VecY,VecZ))*180/pi<<endl;
		
}

void PrimaryGeneratorAction::GetSourceInfo(G4String input)
{
	ifstream ifp2;
	ifp2.open(input);
	string dump1, dump2;
	
	while (!ifp2.eof()) {
	ifp2 >> dump1;
		
	if (dump1 == "pos="){
		ifp2 >> dump2;
		PosX = stod(dump2)*cm;
		ifp2 >> dump2;
		PosY = stod(dump2)*cm;
		ifp2 >> dump2;
		PosZ = stod(dump2)*cm;
	}

	if (dump1 == "vec="){
		ifp2 >> dump2;
		VecX = stod(dump2);
		ifp2 >> dump2;
		VecY = stod(dump2);
		ifp2 >> dump2;
		VecZ = stod(dump2);
	}
	
	if (dump1 == "si1"){
		ifp2 >> dump2;
		ConeAngle = stod(dump2);
		cout<<"ConeAngle: "<<ConeAngle<<endl;
	}
	
	// if (dump1 == "si2"){
	// 	std::vector<G4double> energies;	
	// 	G4double prevCdf(0);
	// 	while(!ifp2.eof()){
	// 	ifp2 >> dump2;
	// 	if (dump2 == "sp2")
	// 	{
	// 		for (int i=0; i<energies.size(); i++){
	// 			ifp2 >> dump2;
	// 			prevCdf += stod(dump2);
	// 			if(i==energies.size()-1) prevCdf = 1.;
	// 			CDF[prevCdf] = energies[i];
	// 		}
	// 		break;
	// 	}
	// 	else 
	// 	{
	// 		energies.push_back(stod(dump2));
	// 		// Eb.push_back(stod(dump2));	
	// 		//cout<<stod(dump2)<<endl;
	// 	}
	// 	}
		
	// }
		
	}
	ifp2.close();

	// for (int i=1; i<Epdf.size(); i++){
	// 	PDF_Up.push_back(pair<G4double, G4double>(Epdf[i],Eb[i]));
	// 	PDF_Dw.push_back(pair<G4double, G4double>(Epdf[i],Eb[i-1]));
	// 	//cout<<PDF_Up[i-1].first<<" "<<PDF_Up[i-1].second<<" "<<PDF_Dw[i-1].second<<endl;
	// }
	
	// sort(PDF_Up.begin(), PDF_Up.end()); 
	// sort(PDF_Dw.begin(), PDF_Dw.end()); 
	
	// G4double cdf=0;
	// for (int i=0; i<PDF_Up.size(); i++){
	// 	cdf=cdf+PDF_Up[i].first;
	// 	CDF.push_back(cdf);

	// 	//cout<<CDF[i]<<" "<<PDF_Up[i].second<<" "<<PDF_Dw[i].second<<endl;
	// }
	
	
	G4ThreeVector V1(0,0,1);
	G4ThreeVector V2(VecX,VecY,VecZ);
	RotAxis=V1.cross(V2);
	RotAngle = V1.angle(V2);
	
}
