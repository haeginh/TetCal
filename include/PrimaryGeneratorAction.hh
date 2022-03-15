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

class TETModelImport;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
	PrimaryGeneratorAction(G4double centerZ);
	virtual ~PrimaryGeneratorAction();

    //GENERAL
  public:
    virtual void   GeneratePrimaries(G4Event* anEvent);

    // G4ParticleGun*  GetParticleGun()            const {return fParticleGun;}
    void SetSpecDir(G4String _specDir) {specDir = _specDir;}
    void SetPeakEnergy(G4int _kVp);
    void SetAngle(G4double _angle) {angle = _angle;}
    void SetRadius(G4double _radius) {radius = _radius;}
    void SetLowerBound(G4double _lb) {lowerBound = _lb-centerZ;}
    void SetUpperBound(G4double _ub) {upperBound = _ub-centerZ;}
    
    G4String GetSourceName() const {
      return std::to_string(kVp)+"kVp_"+std::to_string(int(angle/deg+0.5))+"deg_R"+
             std::to_string(int(radius/cm+0.5))+"cm_H"+
             std::to_string(int(lowerBound/cm+0.5))+"-"+std::to_string(int(upperBound/cm+0.5))+"cm";
      }

    G4double GetFactor() const
    {
      return specSum * (1*m*angle*rad*(upperBound-lowerBound))/cm2;
    }
    
  private:
    G4String             specDir;
    G4int                kVp;
    G4double             angle;
    G4double             radius;
    G4double             lowerBound;
    G4double             upperBound;
    G4double             centerZ;
    G4double             specSum;
    G4ParticleGun*       fParticleGun;
    PrimaryMessenger*    fMessenger;

    std::map<G4double, G4double> cdf;
};

#endif

