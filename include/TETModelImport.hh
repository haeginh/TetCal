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
// TETModelImport.hh
// \file   MRCP_GEANT4/External/include/TETModelImport.hh
// \author Haegin Han
//

#ifndef TETModelImport_h
#define TETModelImport_h 1

#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>

#include "G4UIExecutive.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "G4String.hh"
#include "G4Tet.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Colour.hh"

// *********************************************************************
// This class is to import the phantom data from *.ele, *.node, and
// *.material files.
// -- DataRead: Construct G4Tet by reading data from *.ele and *.node
//              files
// -- MaterialRead: Construct G4Material by reading material data from
//                  *.material file
// -- ColourRead: Construct std::map that contains G4Colour data
//                according to data in colour.dat file
// -- PrintMaterialInformation: Print a table that contains organ ID,
//                              number of tetrahedrons, volume, density,
//                              mass, and name for each organ
// *********************************************************************

class TETModelImport
{
public:
	TETModelImport(G4String phantomName, G4UIExecutive* ui);
    virtual ~TETModelImport() {}

	// get methods
	G4bool        DoseWasOrganized()         { return doseOrganized; }
	std::map<G4int, std::vector<G4int>>
	              GetDoseMap()               { return organ2dose;}
	G4String      GetDoseName(G4int doseID)  { return doseName[doseID];}
	std::map<G4int, G4double> GetDoseMassMap(){ return doseMassMap; }

    G4String      GetPhantomName()           { return phantomName; }
	G4Material*   GetMaterial(G4int idx)     { return materialMap[idx];}
	G4int         GetNumTetrahedron()        { return tetVector.size();}
	G4int         GetMaterialIndex(G4int idx){ return materialVector[idx]; }
	G4Tet*        GetTetrahedron(G4int idx)  { return tetVector[idx]; }
	G4double      GetVolume(G4int idx)       { return volumeMap[idx]; }
	std::map<G4int, G4double> GetMassMap()   { return massMap; }
	std::map<G4int, G4Colour> GetColourMap() { return colourMap; }
	G4ThreeVector GetPhantomSize()           { return phantomSize; }
	G4ThreeVector GetPhantomBoxMin()         { return boundingBox_Min; }
	G4ThreeVector GetPhantomBoxMax()         { return boundingBox_Max; }
	std::map<G4int, G4double> GetRBMratio()  { return rbmRatio;}
	std::map<G4int, G4double> GetBSratio()   { return bsRatio;}
	G4double GetRBMDRF(G4int idx, G4int eIdx){ return rbmDRF[idx][eIdx];}
	G4double GetBSDRF (G4int idx, G4int eIdx){ return bsDRF[idx][eIdx];}
    G4ThreeVector GetAVertex(G4int idx)      { return vertexVector[idx]; }

    std::vector<std::vector<G4int>> GetElements(G4int organID){
        std::vector<std::vector<G4int>> eleVec;
        for(size_t i=0;i<materialVector.size();i++){
            if(materialVector[i]!=organID) continue;
            std::vector<G4int> ele = {eleVector[i][0],
                                      eleVector[i][1],
                                      eleVector[i][2],
                                      eleVector[i][3]};
            std::sort(ele.begin(), ele.end());
            eleVec.push_back(ele);
        }
        return eleVec;
    }

private:

	// private methods
	void DoseRead(G4String);
    void DataRead(G4String, G4String);
	void MaterialRead(G4String);
	void RBMBSRead(G4String);
	void DRFRead(G4String);
	void ColourRead();
	void PrintMaterialInfomation();

	G4int GetID(G4String str) {
		G4String strCut = str.substr(2,str.size()-2);
		size_t pos = 0;
		G4String token;
		while ((pos = strCut.find("_")) != std::string::npos) {
			token = strCut.substr(0, pos);
			break;
		}
		return atoi(token.c_str());
	}

	G4String phantomName;

	G4ThreeVector boundingBox_Min;
	G4ThreeVector boundingBox_Max;
	G4ThreeVector phantomSize;

	std::map<G4int, std::vector<G4int>>   organ2dose;
	std::map<G4int, G4String>  doseName;
	std::map<G4int, G4double>  doseMassMap;
	G4bool                     doseOrganized;

	std::vector<G4ThreeVector> vertexVector;
	std::vector<G4Tet*>        tetVector;
	std::vector<G4int*>        eleVector;
	std::vector<G4int>         materialVector;
	std::map<G4int, G4int>     numTetMap;
	std::map<G4int, G4double>  volumeMap;
	std::map<G4int, G4double>  massMap;
	std::map<G4int, G4Colour>  colourMap;
	std::map<G4int, G4double>  rbmRatio;
	std::map<G4int, G4double>  bsRatio;
	std::map<G4int, std::vector<G4double>> rbmDRF;
	std::map<G4int, std::vector<G4double>> bsDRF;

	std::map<G4int, std::vector<std::pair<G4int, G4double>>> materialIndexMap;
	std::vector<G4int>                                       materialIndex;
	std::map<G4int, G4Material*>                             materialMap;
	std::map<G4int, G4double>                                densityMap;
	std::map<G4int, G4String>                                organNameMap;

};

#endif
