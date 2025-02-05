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

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include <vector>
#include <fstream>
#include <sstream>

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include "G4ParticleGun.hh"

using namespace std;

// class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction(G4String);
   ~PrimaryGeneratorAction();
   void GetSourceInfo(G4String input);
   G4String GetInputName() const {return inp_name;}
   G4double GetEnergy() const {return particleGun->GetParticleEnergy();}

  private:
    void GeneratePrimaries(G4Event* anEvent);
   	void SourceSampling();

  private:
  	G4String inp_name;
    G4ThreeVector direction;
    G4ParticleGun* particleGun;
	
	G4double PosX, PosY, PosZ;
	G4double VecX, VecY, VecZ;
	G4double ConeAngle;
	G4ThreeVector RotAxis;
	G4double RotAngle;
	
	// vector<G4double> Eb;
	// vector<G4double> Epdf;
	// vector<G4double> CDF;
	// map<G4double, G4double> CDF;

	// PrimaryMessenger* fMessenger;
	// G4double coeff;
	
	// vector<pair<G4double, G4double>> PDF_Up;
	// vector<pair<G4double, G4double>> PDF_Dw;
	
};

#endif

