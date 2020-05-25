/*
 * ExternalBeam.cc
 *
 *  Created on: Oct 18, 2019
 *      Author: hhg
 */

#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "../include/ExternalBeam.hh"

ExternalBeam::ExternalBeam(G4ThreeVector _trans)
:beamDir(AP), xHalf(-1), yHalf(-1), zHalf(-1), beamArea(-1), beamCenter(G4ThreeVector()), initChk(false), trans(_trans)
{
    SetDefaultSize();
}

ExternalBeam::~ExternalBeam()
{}

void ExternalBeam::SetBeamDirection(BEAMDIR _dir){
	beamDir = _dir;
    initChk = true;
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
        beamArea = radius*radius*cm2*pi;
		break;
	case ISO:
		beamDirName = "ISO";
        beamArea = radius*radius*cm2*pi;
		break;
	}
}

void ExternalBeam::GetAprimaryPosDir(G4ThreeVector &position, G4ThreeVector &direction)
{
    G4double _r, rand, p1, p2, theta, phi;
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
        _r = radius*sqrt(G4UniformRand())*cm;
		rand = G4UniformRand();
        p1 = _r*cos(rand*2*pi);
        p2 = _r*sin(rand*2*pi);
		theta = G4UniformRand()*2*pi;
		direction  = G4ThreeVector(-1, 0, 0);
		position  = G4ThreeVector(100.*cm, p1, p2);
		direction  = direction.rotateZ(theta);
		position = position.rotateZ(theta);
		break;
	case ISO:
        _r = radius*sqrt(G4UniformRand())*cm;
		rand = G4UniformRand();
        p1 = _r*cos(rand*2*pi);
        p2 = _r*sin(rand*2*pi);
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
    position += beamCenter;
}



