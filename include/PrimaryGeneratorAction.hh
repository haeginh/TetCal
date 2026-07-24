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
#include "PrimaryMessenger.hh"
#include "SourceGenerator.hh"

#include <map>
#include <vector>

class TETModelImport;
class G4ParticleDefinition;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
	PrimaryGeneratorAction(TETModelImport* tetData);
	virtual ~PrimaryGeneratorAction();

    //GENERAL
  public:
    virtual void   GeneratePrimaries(G4Event* anEvent);
    void SetSpectrumSource(G4String inputFile);
    void SetIAEASpectrumSource(G4String inputCommand);

    void SetExternalBeam(){
       if(!fExternal) fExternal = new ExternalBeam();
        fSourceGenerator = fExternal; fSourceGenerator->SetExternal();
    }
    void SetInternalBeam(){
        if(!fInternal) fInternal = new InternalSource(tetData);
        fSourceGenerator = fInternal; fSourceGenerator->SetInternal();
    }
    void SetSurfaceSource(){
        if(!fSurface) fSurface = new SurfaceSource(tetData);
        fSourceGenerator = fSurface; fSourceGenerator->SetInternal();
    }
    void SetRadCodes(G4String inputs){
        RADcodes.clear();
        std::istringstream iss(inputs);
        G4String code;
        while (iss >> code) RADcodes.push_back(code);
    }

    void SetSourceName(G4String _sourceN) {sourceName = _sourceN;}
    G4ParticleGun*  GetParticleGun()            const {return fParticleGun;}
    SourceGenerator* GetSourceGenerator()       const {return fSourceGenerator;}
    ExternalBeam*   GetExternalBeamGenerator()  const {return fExternal;}
    InternalSource* GetInternalBeamGenerator()  const {return fInternal;}
    SurfaceSource*  GetSurfaceSourceGenerator() const {return fSurface;}
    G4String        GetSourceName() const {return sourceName;}
    G4bool          IsSpectrum() const {return spectrumSource;}
    G4bool          IsIAEASpectrum() const {return iaeaSpectrumSource;}

  private:
    struct DiscreteEmission {
        G4ParticleDefinition* particle;
        G4double energy;
        G4double yield;
    };

    struct BetaInterval {
        G4ParticleDefinition* particle;
        G4double energy0;
        G4double energy1;
        G4double weight;
    };

    std::vector<G4String> SplitCommand(G4String input) const;
    G4ParticleDefinition* GetParticleDefinition(G4String particleName) const;
    void LoadIAEADiscrete(G4String inputFile);
    void LoadIAEABeta(G4String inputFile);
    void GenerateIAEADecay(G4Event* anEvent);
    void GenerateOnePrimary(G4Event* anEvent,
                            const G4ThreeVector& position,
                            G4ParticleDefinition* particle,
                            G4double energy);
    G4int SampleMultiplicity(G4double yield) const;
    G4double SampleBetaEnergy(G4ParticleDefinition* particle) const;

    TETModelImport*      tetData;
    G4ParticleGun*       fParticleGun;
    PrimaryMessenger*    fMessenger;
    SourceGenerator*     fSourceGenerator;
    ExternalBeam*        fExternal;
    InternalSource*      fInternal;
    SurfaceSource*       fSurface;
    G4String             sourceName;
    G4bool               spectrumSource;
    G4bool               iaeaSpectrumSource;
    std::map<G4double, G4double> samplingE;
    std::vector<G4String> RADcodes;
    std::vector<DiscreteEmission> iaeaDiscreteEmissions;
    std::vector<BetaInterval>     iaeaBetaIntervals;
    std::map<G4ParticleDefinition*, G4double> iaeaBetaYield;
};

#endif
