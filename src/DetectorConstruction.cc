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

#include "DetectorConstruction.hh"

#include "globals.hh"

#include "G4PhysicalConstants.hh"
#include "G4NistManager.hh"  
#include "G4LogicalVolume.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"

DetectorConstruction::DetectorConstruction(ImportVoxelPhantom* _voxelPhantom, G4String _inp)
:fpWorldPhysical(0) ,voxelPhantom(_voxelPhantom), inp_name(_inp), tableDen(-1.) 
{

    blankAtt = new G4VisAttributes;
    blankAtt->SetVisibility(FALSE);

    fNx = voxelPhantom->GetVoxelResolution(0);
    fNy = voxelPhantom->GetVoxelResolution(1);
    fNz = voxelPhantom->GetVoxelResolution(2);

    fVoxelHalfLengthX = voxelPhantom->GetVoxelSize(0) * 0.5;
    fVoxelHalfLengthY = voxelPhantom->GetVoxelSize(1) * 0.5;
    fVoxelHalfLengthZ = voxelPhantom->GetVoxelSize(2) * 0.5;

    G4cout << fNx << "\t"<< fNy << "\t"<< fNz << "\t"<<G4endl;
    G4cout << fVoxelHalfLengthX << "\t"<< fVoxelHalfLengthY << "\t"<< fVoxelHalfLengthZ << "\t"<<G4endl;
    
    transform = ReadColliAndDAPInfo();
}

DetectorConstruction::~DetectorConstruction()
{}

G4VPhysicalVolume* DetectorConstruction::Construct()
{

    //G4Material* VACUUM = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    G4Material* Pb = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb"); 
    G4Material* carbon;
    if(tableDen>0) carbon = new G4Material("carbon", tableDen,G4NistManager::Instance()->FindOrBuildMaterial("G4_C"));

    // World
    G4VSolid* sol_World = new G4Box("World", 10*m, 10*m, 10*m);
    //G4Material* AIR = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    
    G4Element* elN = new G4Element("Nitrogen","N",7.,14.01*g/mole);
    G4Element* elO = new G4Element("Oxygen","O",8.,16.00*g/mole);
    
    G4Material* AIR = new G4Material("Air",1.30*mg/cm3,2);
    AIR->AddElement(elN, 80.0*perCent);
    AIR->AddElement(elO, 20.0*perCent);
    
    G4LogicalVolume* lv_World = new G4LogicalVolume(sol_World, AIR, "World");
    fpWorldPhysical =
        new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), lv_World, "World", 0, false, -10000);

    // Phantom Geometry
    // 1st - Container for Phantom
    G4String ContainerName("Container");
    G4VSolid* solContainer = new G4Box(ContainerName,
                        fVoxelHalfLengthX * fNx,
                        fVoxelHalfLengthY * fNy,
                        fVoxelHalfLengthZ * fNz);

    G4LogicalVolume* logContainer = new G4LogicalVolume(solContainer,AIR,ContainerName);
    new G4PVPlacement(0, G4ThreeVector(fVoxelHalfLengthX * fNx,fVoxelHalfLengthY * fNy, fVoxelHalfLengthZ * fNz), logContainer, ContainerName, lv_World, false, -10001);

    G4VisAttributes* va_Con = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
    logContainer->SetVisAttributes(va_Con);

    // Replication of Water Phantom Volume.
    // 2nd - Z Slice
    G4String zRepName("RepZ");
    G4VSolid* solZRep = new G4Box(zRepName,fNx*fVoxelHalfLengthX, fNy*fVoxelHalfLengthY, fVoxelHalfLengthZ);
    G4LogicalVolume* logZRep = new G4LogicalVolume(solZRep,AIR,zRepName);
    new G4PVReplica(zRepName,logZRep,logContainer,kZAxis,fNz,fVoxelHalfLengthZ*2.);

    logZRep->SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));

    // 3rd - X Slice
    G4String xRepName("RepX");
    G4VSolid* solXRep = new G4Box(xRepName,fVoxelHalfLengthX, fNy*fVoxelHalfLengthY, fVoxelHalfLengthZ);
    G4LogicalVolume* logXRep = new G4LogicalVolume(solXRep,AIR,xRepName);
    new G4PVReplica(xRepName,logXRep,logZRep,kXAxis,fNx,fVoxelHalfLengthX*2.);


    logXRep->SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));

    // Voxel solid and logical volumes
    // 4th - Y Slice
    G4VSolid* solVoxel = new G4Box("phantom",fVoxelHalfLengthX, fVoxelHalfLengthY,fVoxelHalfLengthZ);
    logicVoxel = new G4LogicalVolume(solVoxel,AIR,"VoxelDetector");


    //
    // Parameterisation for transformation of voxels.
    //  (voxel size is fixed in this example.
    //    e.g. nested parameterisation handles material
    //    and transfomation of voxels.)

    VoxelNestedParameterisation* param =
            new VoxelNestedParameterisation(voxelPhantom);

    new G4PVParameterised("phantom",    // their name
                          logicVoxel, // their logical volume
                          logXRep,      // Mother logical volume
                          kYAxis,       // Are placed along this axis
                          //kUndefined,   // Are placed along this axis
                          fNy,      // Number of cells
                          param);       // Parameterisation.

    //G4cout<<"test: xxx"<<G4endl;
    // Visualization
    //G4VisAttributes* va_World = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    //va_World->SetForceWireframe(true);
    //lv_World->SetVisAttributes(va_World);

    //Collimator
    G4String CollimatorName("Collimator");
    solCollimator = new G4Box(CollimatorName,
                        10.0*cm,
                        0.05*cm,
                        10.0*cm);
        
    G4String BeamPortName("BeamPort");
    solBeamPort = new G4Box(BeamPortName,
                        BeamPortHalfSizeX,
                        0.051*cm,
                        BeamPortHalfSizeZ);
    
    G4String CollimatorName_F("Collimator_F");  
    solCollimator_F = new G4SubtractionSolid(CollimatorName_F,solCollimator, solBeamPort);
                        
    logCollimator = new G4LogicalVolume(solCollimator_F,Pb,CollimatorName_F);
    phyCollimator = new G4PVPlacement(transform, logCollimator, CollimatorName_F, lv_World, false, 10300);

    //cout<<"CollimatorCenter.." << G4ThreeVector(CollimatorCenterX, CollimatorCenterY, CollimatorCenterZ) <<endl;

    G4VisAttributes* va_Col = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
    logCollimator->SetVisAttributes(va_Col);

    //DAP
    // G4String DAPName("DAP");
    // G4VSolid* solDAP = new G4Box(DAPName,
    //                     DAPHalfSizeX,
    //                     DAPHalfSizeY,
    //                     DAPHalfSizeZ);

    // logDAP = new G4LogicalVolume(solDAP,AIR,DAPName);
    // new G4PVPlacement(Rot, G4ThreeVector(DAPCenterX, DAPCenterY, DAPCenterZ), logDAP, DAPName, lv_World, false, 10400);

    //cout<<"CollimatorCenter.." << G4ThreeVector(DAPCenterX, DAPCenterY, DAPCenterZ) <<endl;

    // G4VisAttributes* va_DAP = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
    // logDAP->SetVisAttributes(va_DAP);
    //table
    if(tableDen>0){
        auto table_log = new G4LogicalVolume(new G4Box("table", tableHalfSize.getX(), tableHalfSize.getY(), tableHalfSize.getZ()),
                                             carbon, "table");
        new G4PVPlacement(0, tablePos, table_log, "table", lv_World, false, 40000);
    }
    return fpWorldPhysical;
}

void DetectorConstruction::ConstructSDandField()
{
}

G4Transform3D DetectorConstruction::ReadColliAndDAPInfo()
{
    ifstream ifp2;
    ifp2.open(inp_name);
    string dump1, dump2;
//    tableDen = 1.4*(g/cm3);
//G4double tableThickness = 2*cm;
//tableHalfSize = G4ThreeVector(fVoxelHalfLengthX*fNx,tableThickness*0.5,fVoxelHalfLengthZ*fNz);
//tablePos = tableHalfSize+G4ThreeVector(0, fVoxelHalfLengthY*fNy*2., 0);
    while (!ifp2.eof()) {
    ifp2 >> dump1;
    if (dump1 == "1052") {
        ifp2 >> dump2 >> dump2 >> dump2;
        BeamPortHalfSizeX = stod(dump2)*cm;
        G4cout<<"BeamPortHalfSizeX: "<<BeamPortHalfSizeX<<G4endl;   
    }
    else if (dump1 == "3052") {
        ifp2 >> dump2 >> dump2 >> dump2;
        BeamPortHalfSizeZ = stod(dump2)*cm;
        G4cout<<"BeamPortHalfSizeZ: "<<BeamPortHalfSizeZ<<G4endl;   
    }
    else if (dump1== "4000") {
        ifp2>>dump2>>dump2;
        tableDen = -stod(dump2)*(g/cm3);
        G4cout<<"TableDensity: "<<tableDen/(g/cm3)<<" g/cm3"<<G4endl;
    }
    else if (dump1== "4001") {
        G4double minX, maxX, minY, maxY, minZ, maxZ;
        ifp2>>dump2>>minX>>maxX>>minY>>maxY>>minZ>>maxZ;
        tableHalfSize = 0.5*G4ThreeVector(maxX-minX, maxY-minY, maxZ-minZ)*cm;
        tablePos = 0.5*G4ThreeVector(minX+maxX, minY+maxY, minZ+maxZ)*cm;
        G4cout<<"TableHalfSize(cm): "<<tableHalfSize/cm<<G4endl;
        G4cout<<"TablePos(cm): "<<tablePos/cm<<G4endl;
    }
    // if (dump1 == "1062") {
    //     ifp2 >> dump2;
    //     ifp2 >> dump2;
    //     ifp2 >> dump2;
    //         DAPHalfSizeX = stod(dump2)*cm;
    //         G4cout<<"DAPHalfSizeX: "<<DAPHalfSizeX<<G4endl; 
    // }
    
    // if (dump1 == "2062") {
    //     ifp2 >> dump2;
    //     ifp2 >> dump2;
    //     ifp2 >> dump2;
    //         DAPHalfSizeY = stod(dump2)*cm;
    //         G4cout<<"DAPHalfSizeY: "<<DAPHalfSizeY<<G4endl; 
    // }
    
    // if (dump1 == "3062") {
    //     ifp2 >> dump2;
    //     ifp2 >> dump2;
    //     ifp2 >> dump2;
    //         DAPHalfSizeZ = stod(dump2)*cm;
    //         G4cout<<"DAPHalfSizeZ: "<<DAPHalfSizeZ<<G4endl; 
    // }

    else if (dump1 == "*tr1"){        
        G4ThreeVector trans;
        G4double a, b, c, d, e, f, g, h, i;
        ifp2 >>trans>>a>> b>> c>> d>> e>> f>> g>> h>> i;
        // ifp2.close();
        G4RotationMatrix rot(G4ThreeVector(cos(a*deg), cos(b*deg), cos(c*deg)), 
                    G4ThreeVector(cos(d*deg), cos(e*deg), cos(f*deg)), 
                    G4ThreeVector(cos(g*deg), cos(h*deg), cos(i*deg)));
        
        G4Transform3D transform(rot, trans*cm);
        G4cout<<rot<<G4endl;
        return transform;
    }
    
    // if (dump1 == "*tr2"){
    //     ifp2 >> dump2;
    //     DAPCenterX = stod(dump2)*cm;
    //     ifp2 >> dump2;
    //     DAPCenterY = stod(dump2)*cm;
    //     ifp2 >> dump2;
    //     DAPCenterZ = stod(dump2)*cm;
    // }
    
    else if (dump1 == "vec="){
        ifp2 >> dump2;
        DirectionX = stod(dump2);
        ifp2 >> dump2;
        DirectionY = stod(dump2);
        ifp2 >> dump2;
        DirectionZ = stod(dump2);
    }
    
    }
    ifp2.close();
    
    G4ThreeVector V1 = G4ThreeVector(0,1,0);
    G4ThreeVector V2 = G4ThreeVector(DirectionX,DirectionY,DirectionZ);
    
    G4cout<<"V2: "<< V2 << G4endl;
    G4cout<<"V1: "<< V1 << G4endl;
    
    RotAxis = V2.cross(V1);
    RotAngle = V1.angle(V2);
    Rot = new G4RotationMatrix(RotAxis,RotAngle);
    //cout<<"RotAngle: "<< RotAngle << endl;
    //cout<<"RotatedV: "<< V1.rotate(RotAngle,RotAxis) << endl;
}

void DetectorConstruction::ResetColli(G4String fileName)
{
    delete solCollimator;
    delete solBeamPort;
    delete solCollimator_F;
    delete logCollimator;
    delete phyCollimator;

    inp_name = fileName;
    transform = ReadColliAndDAPInfo();
    //Collimator
    G4String CollimatorName("Collimator");
    solCollimator = new G4Box(CollimatorName,
                        10.0*cm,
                        0.05*cm,
                        10.0*cm);
        
    G4String BeamPortName("BeamPort");
    solBeamPort = new G4Box(BeamPortName,
                        BeamPortHalfSizeX,
                        0.051*cm,
                        BeamPortHalfSizeZ);
    
    G4String CollimatorName_F("Collimator_F");  
    solCollimator_F = new G4SubtractionSolid(CollimatorName_F,solCollimator, solBeamPort);
    G4Material* Pb = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");                    
    logCollimator = new G4LogicalVolume(solCollimator_F,Pb,CollimatorName_F);
    phyCollimator = new G4PVPlacement(transform, logCollimator, CollimatorName_F, fpWorldPhysical->GetLogicalVolume(), false, 10300);
}











