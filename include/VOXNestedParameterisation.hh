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
// $Id: DicomNestedPhantomParameterisation.hh 74857 2013-10-23 07:55:55Z gcosmo $
//
/// \file medical/DICOM/include/DicomNestedPhantomParameterisation.hh
/// \brief Definition of the DicomNestedPhantomParameterisation class
//

#ifndef VOXELNESTEDPARAMETERISATION_HH
#define VOXELNESTEDPARAMETERISATION_HH

#include <vector>
#include <map>

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4VNestedParameterisation.hh"
#include "VOXDetectorConstruction.hh"
#include "VOXModelImport.hh"


class G4VPhysicalVolume;
class G4VTouchable;
class G4VSolid;
class G4Material;
class G4VisAttributes;

// CSG Entities which may be parameterised/replicated
//
class G4Box;
class G4Tubs;
class G4Trd;
class G4Trap;
class G4Cons;
class G4Sphere;
class G4Ellipsoid;
class G4Orb;
class G4Torus;
class G4Para;
class G4Polycone;
class G4Polyhedra;
class G4Hype;

/// Implements a G4VNestedParameterisation

class VOXNestedParameterisation : public G4VNestedParameterisation
{
  public:

    VOXNestedParameterisation(VOXModelImport* );
   ~VOXNestedParameterisation();

    G4Material* ComputeMaterial(G4VPhysicalVolume *currentVol,
                                const G4int repNo,
                                const G4VTouchable *parentTouch );
      // Must cope with parentTouch for navigator's SetupHierarchy

    G4int       GetNumberOfMaterials() const;
    G4Material* GetMaterial(G4int idx) const;
      // Needed to define materials for instances of Nested Parameterisation
      // Current convention: each call should return the materials
      // of all instances with the same mother/ancestor volume

    unsigned int GetMaterialIndex( unsigned int nx, unsigned int ny, unsigned int nz) const;
    unsigned int GetMaterialIndex( unsigned int copyNo) const;

    void ComputeTransformation(const G4int no,
                                     G4VPhysicalVolume *currentPV) const;

    // Additional standard Parameterisation methods,
    // which can be optionally defined, in case solid is used.

    void ComputeDimensions(G4Box &, const G4int,
                                    const G4VPhysicalVolume *) const;


  private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Trd&, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trap&, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Cons&, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Sphere&, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Ellipsoid&, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Orb&, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Torus&, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Para&, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Hype&, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Tubs&, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polycone&, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polyhedra&, const G4int,
                            const G4VPhysicalVolume*) const {}

    using G4VNestedParameterisation::ComputeMaterial;

    void MaterialColourDefinition();

  private:
    VOXModelImport* voxelPhantom;

    G4VisAttributes* blankAtt;
    G4double         fVoxelHalfLengthX, fVoxelHalfLengthY, fVoxelHalfLengthZ;
    G4int			 fNx,fNy,fNz;

};
#endif
