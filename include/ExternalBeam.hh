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
#include "G4Box.hh"
#include "G4LogicalVolumeStore.hh"
#include <vector>

enum BEAMDIR {AP, PA, LLAT, RLAT, ROT, ISO};

class ExternalBeam
{
public:
    ExternalBeam(G4ThreeVector trans);
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
    void SetBeamCenter(G4ThreeVector _center){beamCenter = _center + trans;}
    void SetDefaultSize(){
        G4Box* phantomBox = (G4Box*) G4LogicalVolumeStore::GetInstance()->GetVolume("phantomLogical")->GetSolid();
        xHalf=phantomBox->GetXHalfLength();
        yHalf=phantomBox->GetYHalfLength();
        zHalf=phantomBox->GetZHalfLength();
        radius = zHalf*1.3;
        beamCenter = G4ThreeVector();
        initChk = false;
    }

    G4bool IsInitialized() {return initChk;}

private:
    BEAMDIR  beamDir;
    G4String beamDirName;
    G4double xHalf, yHalf, zHalf;
    G4double radius;
    G4double beamArea;
    G4ThreeVector beamCenter;
    G4bool   initChk;
    G4ThreeVector trans;
};


#endif /* SRC_EXTERNALBEAM_HH_ */
