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

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include <vector>

#include "PrimaryMessenger.hh"
#include "SourceGenerator.hh"

// *********************************************************************
// This is UserPrimaryGeneratorAction, and the source was defined by
// G4GeneralParticleSource class.
// -- GeneratePrimaries: Generate primaries by G4GeneralParticleSource
//                       class.
// *********************************************************************
class TETModelImport;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
	PrimaryGeneratorAction(TETModelImport* tetData);
	virtual ~PrimaryGeneratorAction();

    //GENERAL
  public:
    virtual void   GeneratePrimaries(G4Event* anEvent);
    void           SetExternalBeam()
    	{if(!fExternal) fExternal = new ExternalBeam();
         fSourceGenerator = fExternal; fSourceGenerator->SetExternal();}
    void           SetInternalBeam()
    	{if(!fInternal) fInternal = new InternalSource(tetData);
         fSourceGenerator = fInternal; fSourceGenerator->SetInternal();}
    void SetSourceName(G4String _sourceN) {sourceName = _sourceN;}
    G4ParticleGun*  GetParticleGun()          const {return fParticleGun;}
    SourceGenerator* GetSourceGenerator()      const {return fSourceGenerator;}
    ExternalBeam*   GetExternalBeamGenerator() const {return fExternal;}
    InternalSource* GetInternalBeamGenerator() const {return fInternal;}
    G4String        GetSourceName() const {return sourceName;}

  private:
    TETModelImport*      tetData;
    G4ParticleGun*       fParticleGun;
    PrimaryMessenger* fMessenger;
    SourceGenerator*       fSourceGenerator;
    ExternalBeam*       fExternal;
    InternalSource*       fInternal;
    G4String              sourceName;
};

#endif

