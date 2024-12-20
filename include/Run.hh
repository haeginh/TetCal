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
// TETRun.hh
// \file   MRCP_GEANT4/External/include/TETRun.hh
// \author Haegin Han
//

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4THitsMap.hh"
#include "G4SDManager.hh"
#include "TETModelImport.hh"

typedef std::map<G4int, std::pair<G4double, G4double>> EDEPMAP;

// *********************************************************************
// This is G4Run class that sums up energy deposition from each event.
// The sum of the square of energy deposition was also calculated to
// produce the relative error of the dose.
// -- RecordEvent: Sum up the energy deposition and the square of it.
//                 The sums for each organ were saved as the form of
//                 std::map.
// -- Merge: Merge the data calculated in each thread.
// *********************************************************************

//enum BEAMDIR {AP, PA, LLAT, RLAT, ROT, ISO};

class Run : public G4Run 
{
public:
	Run(TETModelImport* tetData);
	virtual ~Run();

	virtual void RecordEvent(const G4Event*);
    virtual void Merge(const G4Run*);

    EDEPMAP* GetEdepMap()      {return &edepMap;}
    EDEPMAP* GetDosimeterMap() {return &dosimeterMap;}
    std::map<G4int, G4double>* GetSpecMap() {return &specMap;}
    G4String GetParticleName() {return primary;}
    G4String GetBeamDirName()  {return dir;}
    G4double GetBeamEnergy()   {return primaryE;}
    G4double GetBeamArea()     {return beamArea;}
    G4bool   GetIsExternal()   {return isExternal;}


    void SetPrimary(G4String _primary, G4String _dir, G4double _primaryE, G4double _beamArea, G4bool _isExternal)
    {
    	primary = _primary;
    	dir = _dir;
    	primaryE = _primaryE;
    	beamArea = _beamArea;
    	isExternal = _isExternal;
    }

private:
    EDEPMAP edepMap, dosimeterMap;
    std::map<G4int, G4double> specMap;
    G4int   fCollID;
    G4int   fCollID_DRF;
    G4int   fCollID_dosimeter;
    G4String primary;
    G4String dir;
    G4double primaryE;
    G4double beamArea;
    G4bool   isExternal;
    std::map<G4int, std::vector<G4int>>   organ2dose;
	std::map<G4int, G4double>  rbmFactor;
	std::map<G4int, G4double>  bsFactor;
	G4bool doseOrganized;
};

#endif
