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
// TETPrimaryGeneratorAction.hh
// \file   MRCP_GEANT4/External/include/TETPrimaryGeneratorAction.hh
// \author Haegin Han
//

#ifndef TETPrimaryGeneratorAction_h
#define TETPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "TETPrimaryMessenger.hh"
#include "globals.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include <vector>

// *********************************************************************
// This is UserPrimaryGeneratorAction, and the source was defined by
// G4GeneralParticleSource class.
// -- GeneratePrimaries: Generate primaries by G4GeneralParticleSource
//                       class.
// *********************************************************************

enum BEAMDIR {AP, PA, LLAT, RLAT, ROT, ISO};

class TETPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
	TETPrimaryGeneratorAction();
	virtual ~TETPrimaryGeneratorAction();

    //GENERAL
  public:
    virtual void GeneratePrimaries(G4Event* anEvent);
  private:
    G4ParticleGun* fParticleGun;
    TETPrimaryMessenger* fMessenger;
    G4bool   internalSwitch;

    //EXTERNAL
  public:
    void         SetExternalBeam() {internalSwitch = false;}
    void         SetBeamDirection(BEAMDIR _dir);

    G4ParticleGun* GetParticleGun()       const {return fParticleGun;}
    G4String 	   GetBeamDirection() 	  const {return beamDirName;}
    G4double 	   GetbBeamArea() 		  const {return beamArea;}
  private:
    BEAMDIR  beamDir;
    G4String beamDirName;
    G4double xHalf, yHalf, zHalf;
    G4double beamArea;

    //INTERNAL
  public:


  private:

};

#endif

