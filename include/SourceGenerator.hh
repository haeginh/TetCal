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

#include <vector>

class SourceGenerator
{
public:
	virtual ~SourceGenerator() {};
	virtual void GetAprimaryPos(G4ThreeVector &pos)
	{pos = G4ThreeVector(); }
	void SetExternal() {isExternal = true; isInternal = false;}
	void SetInternal() {isInternal = true; isExternal = false;}
	G4bool IsExternal() {return isExternal;}
	G4bool IsInternal() {return isInternal;}
private:
	G4bool isExternal = false;
	G4bool isInternal = false;
};

class    TETModelImport;
typedef  std::pair<G4double, G4int> VOLPICK;
class InternalSource: public SourceGenerator
{
public:
	InternalSource(TETModelImport* tetData);
	virtual ~InternalSource();

	void SetSource(std::vector<G4int> sources);
	void GetAprimaryPos(G4ThreeVector &pos);

	std::vector<G4int> GetSource() 	const {return sourceIDs;}

private:
	G4ThreeVector RandomSamplingInTet(G4Tet* tet);

private:
    std::vector<G4int>    sourceIDs;
    TETModelImport*       tetData;
    std::vector<VOLPICK>  tetPick;
};



#endif /* SRC_EXTERNALBEAM_HH_ */
