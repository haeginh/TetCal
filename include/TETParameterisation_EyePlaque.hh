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
// TETParameterisation_EyePlaque.hh
// \file   MRCP_GEANT4/External/include/TETParameterisation_EyePlaque.hh
// \author Haegin Han
//

#ifndef TETParameterisation_EyePlaque_h
#define TETParameterisation_EyePlaque_h 1


#include "globals.hh"
#include "G4VPVParameterisation.hh"
#include "G4VSolid.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4Tet.hh"

#include <map>

class G4VPhysicalVolume;

// *********************************************************************
// This class defines the phantom geometry by using G4PVParameterisation
// class.
// -- ComputeSolid: return the G4Tet* for each element
// -- ComputeMaterial: return the G4Material* corresponding to each organ,
//                     and set the colours for visualization purposes
// *********************************************************************

class TETParameterisation_EyePlaque : public G4VPVParameterisation
{
  public:
    TETParameterisation_EyePlaque(std::vector<G4Tet*>, std::vector<G4int>, std::map<G4int, G4Material*>);
    virtual ~TETParameterisation_EyePlaque();
    
    virtual G4VSolid* ComputeSolid(
    		       const G4int copyNo, G4VPhysicalVolume* );
    
    virtual void ComputeTransformation(
                   const G4int,G4VPhysicalVolume*) const;

    virtual G4Material* ComputeMaterial(const G4int copyNo,
                                        G4VPhysicalVolume* phy,
                                        const G4VTouchable*);

  private:
    std::vector<G4Tet*>                tetVec;
    std::vector<G4int>                 idVec;
    std::map<G4int, G4Material*>       matMap;
    std::map<G4int, G4VisAttributes*>  visAttMap;
};

#endif
