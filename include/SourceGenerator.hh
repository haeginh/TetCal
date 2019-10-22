/*
 * ExternalBeam.hh
 *
 *  Created on: Oct 18, 2019
 *      Author: hhg
 */

#ifndef SRC_BEAMGENERATOR_HH_
#define SRC_BEAMGENERATOR_HH_

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Tet.hh"
#include "G4UTet.hh"

#include <vector>

class SourceGenerator
{
public:
	virtual ~SourceGenerator() {};
	virtual void GetAprimaryPosDir(G4ThreeVector &pos, G4ThreeVector &dir)
	{pos = G4ThreeVector(); dir = G4ThreeVector(0., 0., -1.);}
	void SetExternal() {isExternal = true; isInternal = false;}
	void SetInternal() {isInternal = true; isExternal = false;}
	G4bool IsExternal() {return isExternal;}
	G4bool IsInternal() {return isInternal;}
private:
	G4bool isExternal = false;
	G4bool isInternal = false;
};

enum BEAMDIR {AP, PA, LLAT, RLAT, ROT, ISO};

class ExternalBeam: public SourceGenerator
{
public:
	ExternalBeam();
	virtual ~ExternalBeam();

	void     SetBeamDirection(BEAMDIR _dir);
	void GetAprimaryPosDir(G4ThreeVector &pos, G4ThreeVector &dir);

	G4String GetBeamDirection() 	  const {return beamDirName;}
    G4double GetBeamArea() 	  	      const {return beamArea;}

private:
    BEAMDIR  beamDir;
    G4String beamDirName;
    G4double xHalf, yHalf, zHalf;
    G4double beamArea;
};

class    TETModelImport;
typedef  std::pair<G4double, G4int> VOLPICK;
class InternalSource: public SourceGenerator
{
public:
	InternalSource(TETModelImport* tetData);
	virtual ~InternalSource();

	void SetSource(std::vector<G4int> sources);
	void GetAprimaryPosDir(G4ThreeVector &pos, G4ThreeVector &dir);

	std::vector<G4int> GetSource() 	const {return sourceIDs;}

private:
	G4ThreeVector RandomSamplingInTet(G4Tet* tet);

private:
    std::vector<G4int>    sourceIDs;
    TETModelImport*       tetData;
    std::vector<VOLPICK>  tetPick;
};

/*class SurfaceSource: public SourceGenerator
{
public:
	SurfaceSource(TETModelImport* tetData);
	virtual ~SurfaceSource();

	void SetSource(std::vector<G4int> sources);
	void GetAprimaryPosDir(G4ThreeVector &pos, G4ThreeVector &dir);

	std::vector<G4int> GetSource() 	const {return sourceIDs;}

private:
	G4ThreeVector RandomSamplingInTriangle(G4Tet* tet);

private:
    std::vector<G4int>    sourceIDs;
    TETModelImport*       tetData;
    std::vector<VOLPICK>  tetPick;
};*/


#endif /* SRC_EXTERNALBEAM_HH_ */
