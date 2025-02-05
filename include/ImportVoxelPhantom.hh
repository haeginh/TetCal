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

#ifndef INCLUDE_IMPORTVOXELPHANTOM_HH_
#define INCLUDE_IMPORTVOXELPHANTOM_HH_

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>
#include <cstring>

#include "G4ios.hh"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Element.hh"

using namespace std;

class ImportVoxelPhantom {
public:
	ImportVoxelPhantom(G4String);
	virtual ~ImportVoxelPhantom();


	void ImportPhantomInput();
	void ImportPhantomVoxelData();
	void ImportPhantomVoxelMaterial();
	void ImportPhantomVoxelVolume();

	G4String GetInputName (){return inp_name;}
	G4String GetLatticeName (){return lat_name;}
	G4int GetVoxelResolution(G4int idx) {return voxelResolution[idx];}
	G4double GetVoxelSize(G4int idx) {return voxelSize[idx];}
	G4Material* GetVoxelMaterial(G4int idx) {return materialMap[idx];}
	G4double GetOrganVolume(G4int idx) {return organVolume[idx];}
	G4double GetOrganMass(G4int idx) {return organVolume[idx] * Density[idx] * g/cm3;}
	G4double GetDAPMass() {return DAPHalfSizeX*DAPHalfSizeY*DAPHalfSizeZ*8*0.0013* g/cm3;}
	G4int GetMaterialIndex(G4int idx) {return materialIndex[idx];}
	G4int GetVoxelMaterialSize() {return (G4int) materialMap.size();}
	G4int GetVoxelData(G4int idx, G4int idy, G4int idz) {return voxelData[idx][idy][idz];}
	std::map<G4int, G4String> GetOrganNameAll() {return OrganNameAll;}
	G4int GetNumVoxel(G4int id){return numVoxel[id];}


private:
	std::map<G4int, G4double> organVolume;
	std::map<G4int, G4int> numVoxel;
	std::map<G4int, std::vector<std::pair<G4int, G4double> > > materialIndexMap;
	std::vector<G4int> materialIndex;
	std::map<G4int, G4Material* > materialMap;
	std::map<G4int, G4double> densityMap;
	std::map<G4int, G4String> MatNameMap;
	std::map<G4int, G4String> OrganNameAll;


	std::map<G4int, G4double> RBM;
	std::map<G4int, G4double> BS;

	std::vector<G4int> RBMIndex;

	G4ThreeVector voxelSize;
	std::vector<G4int> voxelResolution;
	G4ThreeVector phantomSize;

	G4int*** voxelData;
	G4int VoxIDtoMatID[2][256];
	G4int CompID[256]; 
	G4double Density[256];
	G4double DAPHalfSizeX;
	G4double DAPHalfSizeY;
	G4double DAPHalfSizeZ;
	G4String inp_name;
	G4String lat_name;
	G4int Gender; //male:0, female:1

};

#endif /* INCLUDE_IMPORTVOXELPHANTOM_HH_ */
