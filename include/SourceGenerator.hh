
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

#include <vector>

class SourceGenerator
{
public:
    virtual ~SourceGenerator() {}
    virtual void GetAprimaryPos(G4ThreeVector &pos)
    {pos = G4ThreeVector();}
    void SetExternal() {isExternal = true; isInternal = false;}
    void SetInternal() {isInternal = true; isExternal = false;}
    G4bool IsExternal() {return isExternal;}
    G4bool IsInternal() {return isInternal;}
private:
    G4bool isExternal = false;
    G4bool isInternal = false;
};

enum BEAMDIR {AP, PA, LLAT, RLAT, ROT, ISO};

class    VOXModelImport;
typedef  std::tuple<G4int, G4int, G4int> VOX;
class InternalSource: public SourceGenerator
{
public:
    InternalSource(VOXModelImport* voxData);
    virtual ~InternalSource();

    void SetSource(std::vector<G4int> sources);
    void GetAprimaryPos(G4ThreeVector &pos);

    std::vector<G4int> GetSource() 	const {return sourceIDs;}

private:
    G4ThreeVector RandomSamplingInAVoxel(VOX vox);

private:
    std::vector<G4int>    sourceIDs;
    VOXModelImport*       voxData;
    std::map<G4int,std::vector<VOX>>  voxPick;
    G4ThreeVector         base;
    G4bool                isRBM;
    std::vector<std::pair<G4double, G4int>> spPick;
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
