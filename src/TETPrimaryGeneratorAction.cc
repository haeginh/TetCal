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


#include "TETPrimaryGeneratorAction.hh"
#include "G4LogicalVolumeStore.hh"
#include <fstream>
#include "G4Box.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"

TETPrimaryGeneratorAction::TETPrimaryGeneratorAction()
:G4VUserPrimaryGeneratorAction(), fParticleGun(0), beamDir(AP), xHalf(-1), yHalf(-1), zHalf(-1), beamArea(-1),
 internalSwitch(false)
{
	fParticleGun = new G4ParticleGun(1);

	//initialization for phantom box
	G4Box* phantomBox = (G4Box*) G4LogicalVolumeStore::GetInstance()->GetVolume("phantomLogical")->GetSolid();
	xHalf=phantomBox->GetXHalfLength();
	yHalf=phantomBox->GetYHalfLength();
	zHalf=phantomBox->GetZHalfLength();
	fMessenger = new TETPrimaryMessenger(this);
}

TETPrimaryGeneratorAction::~TETPrimaryGeneratorAction()
{
	delete fParticleGun;
}

void TETPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	G4ThreeVector direction, position;
	G4double radius, rand, p1, p2, theta, phi;
	switch(beamDir)
	{
	case AP:
		direction  = G4ThreeVector(0, 1, 0);
		position.setX(-xHalf+2*xHalf*G4UniformRand());
		position.setY(-200.*cm);
		position.setZ(-zHalf+2*zHalf*G4UniformRand());
		break;
	case PA:
		direction  = G4ThreeVector(0, -1, 0);
		position.setX(-xHalf+2*xHalf*G4UniformRand());
		position.setY(200.*cm);
		position.setZ(-zHalf+2*zHalf*G4UniformRand());
		break;
	case LLAT:
		direction  = G4ThreeVector(-1, 0, 0);
		position.setX(200*cm);
		position.setY(-yHalf+2*yHalf*G4UniformRand());
		position.setZ(-zHalf+2*zHalf*G4UniformRand());
		break;
	case RLAT:
		direction  = G4ThreeVector(1, 0, 0);
		position.setX(-200*cm);
		position.setY(-yHalf+2*yHalf*G4UniformRand());
		position.setZ(-zHalf+2*zHalf*G4UniformRand());
		break;
	case ROT:
		radius = 100*sqrt(G4UniformRand())*cm;
		rand = G4UniformRand();
		p1 = radius*cos(rand*2*pi);
		p2 = radius*sin(rand*2*pi);
		theta = G4UniformRand()*2*pi;
		direction  = G4ThreeVector(-1, 0, 0);
		position  = G4ThreeVector(100.*cm, p1, p2);
		direction  = direction.rotateZ(theta);
		position = position.rotateZ(theta);
		break;
	case ISO:
		radius = 100*sqrt(G4UniformRand())*cm;
		rand = G4UniformRand();
		p1 = radius*cos(rand*2*pi);
		p2 = radius*sin(rand*2*pi);
		theta = G4UniformRand()*2*pi;
		phi = acos(G4UniformRand()*2.-1.);
		direction = G4ThreeVector(0, 0, -1.);
		position = G4ThreeVector(p1, p2, 200.*cm);
		direction = direction.rotateY(phi);
		position = position.rotateY(phi);
		direction = direction.rotateZ(theta);
		position = position.rotateZ(theta);
		break;
	}


	fParticleGun->SetParticlePosition(position);
	fParticleGun->SetParticleMomentumDirection(direction);
	fParticleGun->GeneratePrimaryVertex(anEvent);
}

void TETPrimaryGeneratorAction::SetBeamDirection(BEAMDIR _dir){
	beamDir = _dir;
	switch(beamDir){
	case AP:
		beamDirName = "AP";
		beamArea = xHalf*zHalf*4.;
		break;
	case PA:
		beamDirName = "PA";
		beamArea = xHalf*zHalf*4.;
		break;
	case LLAT:
		beamDirName = "LLAT";
		beamArea = yHalf*zHalf*4.;
		break;
	case RLAT:
		beamDirName = "RLAT";
		beamArea = yHalf*zHalf*4.;
		break;
	case ROT:
		beamDirName = "ROT";
		beamArea = 10000*cm2*pi;
		break;
	case ISO:
		beamDirName = "ISO";
		beamArea = 10000*cm2*pi;
		break;
	}
}

