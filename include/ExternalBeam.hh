/*
 * ExternalBeam.hh
 *
 *  Created on: Oct 18, 2019
 *      Author: hhg
 */

#ifndef SRC_EXTERNALBEAM_HH_
#define SRC_EXTERNALBEAM_HH_

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Tet.hh"

#include <vector>

enum BEAMDIR {AP, PA, LLAT, RLAT, ROT, ISO};

class ExternalBeam
{
public:
	ExternalBeam();
	virtual ~ExternalBeam();

	void     SetBeamDirection(BEAMDIR _dir);
	void GetAprimaryPosDir(G4ThreeVector &pos, G4ThreeVector &dir);

	G4String GetBeamDirection() 	  const {return beamDirName;}
    G4double GetBeamArea() 	  	      const {return beamArea;}

    //For beam area setting
    void SetBeamXYZ(G4ThreeVector size){
        xHalf = size.getX()*0.5;
        yHalf = size.getY()*0.5;
        zHalf = size.getZ()*0.5;
        initChk = false;
    }
    void SetBeamRadius(G4double _r) {radius = _r; initChk = false;}
    void SetBeamCenter(G4ThreeVector _center){beamCenter = _center;}

    G4bool IsInitialized() {return initChk;}

private:
    BEAMDIR  beamDir;
    G4String beamDirName;
    G4double xHalf, yHalf, zHalf;
    G4double radius;
    G4double beamArea;
    G4ThreeVector beamCenter;
    G4bool   initChk;
};


#endif /* SRC_EXTERNALBEAM_HH_ */
