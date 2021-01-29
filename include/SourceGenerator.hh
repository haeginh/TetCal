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
#include "G4RandomDirection.hh"
#include <vector>

class SourceGenerator
{
public:
    virtual ~SourceGenerator() {}
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
    void SetSource(std::vector<G4int> sources, std::vector<G4double> fractions);
	void GetAprimaryPosDir(G4ThreeVector &pos, G4ThreeVector &dir);

	std::vector<G4int> GetSource() 	const {return sourceIDs;}

private:
	G4ThreeVector RandomSamplingInTet(G4Tet* tet);

private:
    std::vector<G4int>    sourceIDs;
    TETModelImport*       tetData;
    std::vector<VOLPICK>  tetPick;
};

class LungSource: public SourceGenerator
{
public:
    LungSource(){}
    virtual ~LungSource(){}

    void SetSource(G4String file){
        sources.clear();
        std::ifstream ifs(file);
        if(!ifs.is_open()){
            G4cerr<<file<<" is not open"<<G4endl;
            exit(300);
        }
        G4double xPos,yPos,zPos;
        for(G4int i=0;i<100000000;i++){
            ifs>>xPos>>yPos>>zPos;
            sources.push_back(G4ThreeVector(xPos, yPos, zPos)*cm);
        }
    }
    void GetAprimaryPosDir(G4ThreeVector &pos, G4ThreeVector &dir){
        pos = sources[floor(G4UniformRand()*100000000)];
        dir = G4RandomDirection();
    }

private:
    std::vector<G4ThreeVector> sources;
};

typedef std::pair<G4double, std::tuple<G4int, G4int, G4int>> TRIPICK;
class SurfaceSource: public SourceGenerator
{
public:
	SurfaceSource(TETModelImport* tetData);
	virtual ~SurfaceSource();

	void SetSource(std::vector<G4int> sources);
	void GetAprimaryPosDir(G4ThreeVector &pos, G4ThreeVector &dir);

	std::vector<G4int> GetSource() 	const {return sourceIDs;}

private:
    G4double      CalculateTriangleArea(G4ThreeVector a, G4ThreeVector b, G4ThreeVector c);
    G4ThreeVector RandomSamplingInTriangle(std::tuple<G4int, G4int, G4int>);

private:
    std::vector<G4int>    sourceIDs;
    TETModelImport*       tetData;
    std::vector<TRIPICK>  facePick;
};


#endif /* SRC_EXTERNALBEAM_HH_ */
