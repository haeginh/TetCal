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


#ifndef DRFScorer_h
#define DRFScorer_h 1

#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"

#include "VOXModelImport.hh"

class DRFScorer : public G4VPrimitiveScorer
{
   public: // with description
      DRFScorer(G4String name,VOXModelImport*);
      virtual ~DRFScorer();

  protected: // with description
      virtual G4int GetIndex(G4Step*);
      virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

  public:
      virtual void Initialize(G4HCofThisEvent*);
      virtual void EndOfEvent(G4HCofThisEvent*);
      virtual void clear();

  private:

      VOXModelImport* PhantomData;
      G4int FindIndexfromEnergyBin(G4double energy);
      G4double GetRBMdose(G4double energy, G4double cellFlux, G4int organID);
      G4double GetBSdose(G4double energy, G4double cellFlux, G4int organID);

      G4int HCID;
      G4THitsMap<G4double>* EvtMap;
      std::vector<G4double> energyBin;
      std::map<G4int, G4double> RBMratio;
      std::map<G4int, G4double> BSratio;

      G4ParticleDefinition* gamma;
};
#endif

