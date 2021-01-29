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
// The code was written by : hhg
//
//
#include "LungParallelDetCon.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "Randomize.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4UIcommand.hh"


#include <cmath>

G4ThreadLocal G4bool LungParallelDetCon::fSDConstructed = false;

LungParallelDetCon
::LungParallelDetCon(G4String parallelWorldName, TETModelImport* _tetData)
:G4VUserParallelWorld(parallelWorldName),fConstructed(false),tetData(_tetData),
 bBoxLogical(0), volChkName("")
{
    fMessenger = new LungDetMessenger(this);
    bbox.first = G4ThreeVector(1, 1, 1);
	bbox.second = G4ThreeVector(-1, -1, -1);
	box_center = G4ThreeVector();

	for(G4int i=0;i<16;i++){
		for(G4int j=0;j<10;j++){
			bb_diam[i][j]=0;
		}
	}
}

LungParallelDetCon::~LungParallelDetCon()
{
    delete fMessenger;
}

void LungParallelDetCon::Construct()
{
	if(fConstructed) return;
	fConstructed = true;

	//
	// World
	//
	G4VPhysicalVolume* ghostWorld = GetWorld();
	G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();

	//
	// material defined in the mass world
	//
	for(int i=160;i<169;i++){ //bb
        mat_map[i] = tetData->GetMaterial(800+i-160);
	}
    mat_map[146] = tetData->GetMaterial(14000); //air

	//
	// import lung data (info file / diameters)
	//
	ImportInfoData();
	SetBBdiameters();
	PrintLungData();

	//
	// set the bounding box for the computation time (margin: 5 cm)
	//
	G4ThreeVector box_dim = bbox.first - bbox.second;
	box_center = (bbox.first + bbox.second) * 0.5;

	G4VSolid* bBox = new G4Box("bBox",box_dim.x() * .5 + 3. *cm,
			box_dim.y() * .5 + 3. *cm,
			box_dim.z() * .5 + 3. *cm);

	bBoxLogical = new G4LogicalVolume(bBox,0,"bBox");
	new G4PVPlacement(0,box_center,bBoxLogical,"bBox",worldLogical,false,2000);

	//
	// construct the BBs
	//
	G4int pre_mid_gen(0);
	G4int num_15gen = pow(2, 15)-1;
	for(G4int mid_id=2; mid_id<num_15gen;mid_id++){
		G4int mid_gen = (G4int)log2(mid_id);
		if(pre_mid_gen!=mid_gen)
			G4cout<<"gen-"<<mid_gen+1<<"..............";

		std::vector<BB> fBB_vec;

		fBB_vec.push_back(bb_vec[(G4int)(mid_id*0.5)]);
		fBB_vec.push_back(bb_vec[mid_id]);
		fBB_vec.push_back(bb_vec[mid_id*2]);
		fBB_vec.push_back(bb_vec[mid_id*2+1]);

		ConstructBBUnit(fBB_vec, mid_gen);

		if(pre_mid_gen!=mid_gen)
			G4cout<<"constructed"<<G4endl;

		pre_mid_gen=mid_gen;
	}

	//
	//draw bb
	//
    std::map<int, G4double> vol_data;

    G4cout<< phys_vec.size() << " physical volumes were constructed.." << G4endl;

    if(!volChkName.empty()){
        G4int count(0);
        for(auto phyV:phys_vec){
            G4double vol = phyV->GetLogicalVolume()->GetSolid()->GetCubicVolume();
            G4cout<<++count<<"	"<<phyV->GetCopyNo()
                               <<"	"<<phyV->GetLogicalVolume()->GetName()
                               <<"	"<<phyV->GetLogicalVolume()->GetMaterial()->GetName()
                               <<"	"<<vol/mm3<<" mm3"<<G4endl;
            vol_data[phyV->GetCopyNo()]+=vol;
        }

        std::ofstream ofs_vol;
        ofs_vol.open(G4String(volChkName));
        G4double prevVol(0.);
        for(auto data:vol_data){
            if(data.first%1000==146) prevVol=0;
            ofs_vol<<data.first<<": "<<data.second-prevVol<<G4endl;
            prevVol = data.second;
        }
        ofs_vol<<"BB-bas(1163+1164): "<<vol_data[1164]-vol_data[1162] <<" mm3"<<G4endl;
        ofs_vol<<"BB-sec(1164+1165): "<<vol_data[1165]-vol_data[1163] <<" mm3"<<G4endl;
        ofs_vol<<"bb-sec(2165)     : "<<vol_data[2165]-vol_data[2164] <<" mm3"<<G4endl;
        ofs_vol.close();
	}

}
void LungParallelDetCon::ImportInfoData()
{
	//read file
	std::ifstream ifs;
    ifs.open((tetData->GetPhantomName()+".lung").c_str());
    if(!ifs) G4cerr<<tetData->GetPhantomName()<<".lung file is not open"<<G4endl;

	G4int dump;
	G4bool chk;
	G4ThreeVector bb_center;

	ifs>>dump>>bb_center>>chk;
	bb_center = bb_center * cm;
	bb_vec.push_back(BB(bb_center,chk));
	bbox.first = bb_center;
	bbox.second = bb_center;

	while(!ifs.eof()){
		ifs>>dump>>bb_center>>chk;
		bb_center = bb_center * cm;

		if(bb_center.x() > bbox.first.x()) bbox.first.setX(bb_center.x());
		if(bb_center.y() > bbox.first.y()) bbox.first.setY(bb_center.y());
		if(bb_center.z() > bbox.first.z()) bbox.first.setZ(bb_center.z());

		if(bb_center.x() < bbox.second.x()) bbox.second.setX(bb_center.x());
		if(bb_center.y() < bbox.second.y()) bbox.second.setY(bb_center.y());
		if(bb_center.z() < bbox.second.z()) bbox.second.setZ(bb_center.z());

		bb_vec.push_back(BB(bb_center,chk));
	}
}

void LungParallelDetCon::SetBBdiameters()
{
	//airway diameter data
    G4String diamFile = tetData->GetPhantomName() + ".lungDiam";
	std::ifstream ifs(diamFile.c_str());
    if(!ifs) G4cerr<<tetData->GetPhantomName()<<".lung file is not open"<<G4endl;

	G4double douTemp;
	G4int i(0);

	while(ifs>>douTemp){
		if((G4int)(i/10)<8)
			bb_diam[(G4int)(i/10)][i%10] = douTemp * m;
		else
			bb_diam[(G4int)((i-80)/8) + 8][i%8] = douTemp * m;
		i++;
	}
	for(int j=0;j<8;j++){
		bb_diam[16][j] = bb_diam[15][j];
	}
}

void LungParallelDetCon::PrintLungData()
{
	G4int num_in_gen[15];

	for(G4int gen=0;gen<15;gen++){
		num_in_gen[gen] = 0;
		for(G4int i=pow(2, gen+1);i<pow(2, gen+2);i++){
			if(bb_vec[i].second)
				num_in_gen[gen]++;
		}
	}

	G4cout<<G4endl<<"*****************************************Lung Data*******************************************"<<G4endl;
	G4cout<<"[mm]	#";
	for(G4int j = 0;j<10;j++)
		G4cout<<"	layer"<<j;
	G4cout<<G4endl;

	for(int i=0;i<15;i++){
		G4cout<<"gen-"<<i+1<<"	"<<num_in_gen[i]<<"	";
		for(G4int j=0;j<10;j++) G4cout<<bb_diam[i][j]<<"	";
		G4cout<<G4endl;
	}
	G4cout<<"gen-16"<<"		";
	for(G4int j=0;j<10;j++) G4cout<<bb_diam[15][j]<<"	";
	G4cout<<G4endl<<"*********************************************************************************************"<<G4endl;
}

void LungParallelDetCon::ConstructBBUnit(std::vector<BB> fBB_vec, G4int mid_gen)
{
	if(!fBB_vec[0].second) return ; //if there is no grandmother
	if(!fBB_vec[1].second) return ; //it there is no mother
	if((!fBB_vec[2].second) && (!fBB_vec[3].second)) return;

	//Set the number of layers
	G4int num_Layers(10);
	G4int BB1000_bb2000(1000);
	if(mid_gen > 7){
		num_Layers = 8;
		BB1000_bb2000 = 2000;
	}

	//center spheres for each layer (solid)
	std::vector<G4VSolid*> center_vec;
	for(G4int l=0;l<num_Layers;l++)
		center_vec.push_back(new G4Orb("center", bb_diam[mid_gen][l]*0.5));

	//daughter branches (solid)
	std::vector<G4VSolid*> dBranches; //dir1_layer0, dir1_layer1, dir1_layer2.....dir2_layer0, dir2_layer1...
	std::vector<G4double> h2_vec;
	std::vector<G4double> theta2_vec;

	for(G4int gen=2;gen<4;gen++){
		if(!fBB_vec[gen].second) continue;

		//G4double d(fBB_vec[1].first.howNear(fBB_vec[gen].first));
		G4double d= sqrt(pow(fBB_vec[1].first.x()-fBB_vec[gen].first.x(), 2)+
						 pow(fBB_vec[1].first.y()-fBB_vec[gen].first.y(), 2)+
						 pow(fBB_vec[1].first.z()-fBB_vec[gen].first.z(), 2));

		G4double box_hl = (d>bb_diam[mid_gen][0])? d: bb_diam[mid_gen][0];
		G4VSolid* box = new G4Box("box to cut", box_hl, box_hl, box_hl);

		for(G4int l=0;l<num_Layers;l++){
			G4double theta = asin((bb_diam[mid_gen][l]-bb_diam[mid_gen+1][l])*.5/d);
			theta2_vec.push_back(theta);
			G4double r_l = bb_diam[mid_gen][l]*.5*cos(theta);
			G4double r_u = bb_diam[mid_gen+1][l]*.5*cos(theta);
			G4double h = d*cos(theta)*cos(theta);
			h2_vec.push_back(h);
			G4VSolid* con = new G4Cons("dau_cone_temp", 0., r_l, 0., r_u, h*0.5, 0, 2*pi);

			dBranches.push_back(new G4SubtractionSolid ("daughter_cone", con, box,
								0, G4ThreeVector(0., 0., d+box_hl-bb_diam[mid_gen+1][0]*0.5-0.5*h)));
		}
	}

	//mother branches (solid)
	std::vector<G4VSolid*> mBranches;
	std::vector<G4double> theta_vec;
	std::vector<G4double> h_vec;

	for(G4int i=0;i<1;i++){
		G4double d = sqrt(pow(fBB_vec[1].first.x()-fBB_vec[0].first.x(), 2)+
						  pow(fBB_vec[1].first.y()-fBB_vec[0].first.y(), 2)+
						  pow(fBB_vec[1].first.z()-fBB_vec[0].first.z(), 2));

		G4double box_hl = (d>bb_diam[mid_gen-1][0])? d:bb_diam[mid_gen-1][0];
		G4VSolid* box = new G4Box("box to cut", box_hl, box_hl, box_hl);

		for(G4int l=0;l<num_Layers;l++){
			G4double theta = asin((bb_diam[mid_gen-1][l]-bb_diam[mid_gen][l])*.5/d);
			theta_vec.push_back(theta);
			G4double r_u = bb_diam[mid_gen-1][l]*.5*cos(theta);
			G4double r_l = bb_diam[mid_gen][l]*.5*cos(theta);
			G4double h = d*cos(theta)*cos(theta);
			h_vec.push_back(h);
			G4VSolid* con = new G4Cons("mom_cone_temp", 0., r_l, 0., r_u, h*0.5, 0, 2*pi);

			mBranches.push_back(new G4SubtractionSolid ("mother_cone", con, box,
								0,G4ThreeVector(0., 0.,box_hl+bb_diam[mid_gen][0]*0.5-0.5*h)));
		}
	}

	//generate union geometry
	G4ThreeVector cone_vec1 = fBB_vec[0].first - fBB_vec[1].first;
	G4ThreeVector cone_unit_vec1 = cone_vec1 / cone_vec1.r();
	G4ThreeVector unitVectorZ(0., 0., 1.);
	G4double rotation_angle1 = acos(unitVectorZ.dot(cone_unit_vec1));
	G4ThreeVector rotation_axis1 = cone_unit_vec1.cross(unitVectorZ);
	G4RotationMatrix* rot1 = new G4RotationMatrix;
	rot1->rotate(rotation_angle1, rotation_axis1);

	G4ThreeVector cone_vec2 = fBB_vec[2].first - fBB_vec[1].first;
	G4ThreeVector cone_unit_vec2 = cone_vec2 / cone_vec2.r();
	G4double rotation_angle2 = acos(unitVectorZ.dot(cone_unit_vec2));
	G4ThreeVector rotation_axis2 = cone_unit_vec2.cross(unitVectorZ);
	G4RotationMatrix* rot2 = new G4RotationMatrix;
	rot2->rotate(rotation_angle2, rotation_axis2);

	G4ThreeVector cone_vec3 = fBB_vec[3].first - fBB_vec[1].first;
	G4ThreeVector cone_unit_vec3 = cone_vec3 / cone_vec3.r();
	G4double rotation_angle3 = acos(unitVectorZ.dot(cone_unit_vec3));
	G4ThreeVector rotation_axis3 = cone_unit_vec3.cross(unitVectorZ);
	G4RotationMatrix* rot3 = new G4RotationMatrix;
	rot3->rotate(rotation_angle3, rotation_axis3);

	G4LogicalVolume* pre_unit_log=0;
	for(G4int l=0;l<num_Layers;l++){
		//construct union solid
		G4double moving_dist1 = h_vec[l]*0.5-bb_diam[mid_gen][l]*0.5*sin(theta_vec[l]);
		G4VSolid* unit = new G4UnionSolid("center+mom",
										  center_vec[l],
										  mBranches[l],
										  rot1,
										  cone_unit_vec1*moving_dist1);
		if(fBB_vec[2].second && fBB_vec[3].second){ //--if there are both of the daughters
			G4double moving_dist2 = h2_vec[l]*0.5+bb_diam[mid_gen][l]*0.5*sin(theta2_vec[l]);
			G4double moving_dist3 = h2_vec[l+num_Layers]*0.5+bb_diam[mid_gen][l]*0.5*sin(theta2_vec[l+num_Layers]);
			unit = new G4UnionSolid("center+mom+1dau_temp", unit, dBranches[l], rot2, cone_unit_vec2*moving_dist2);
			unit = new G4UnionSolid("center+mom+2dau", unit, dBranches[l+num_Layers], rot3,  cone_unit_vec3*moving_dist3);
		}
		else if(fBB_vec[2].second){
			G4double moving_dist2 = h2_vec[l]*0.5+bb_diam[mid_gen][l]*0.5*sin(theta2_vec[l]);
			unit = new G4UnionSolid("center+mom+dau1", unit, dBranches[l], rot2,cone_unit_vec2*moving_dist2);
		}
		else if(fBB_vec[3].second){
			G4double moving_dist3 = h2_vec[l]*0.5+bb_diam[mid_gen][l]*0.5*sin(theta2_vec[l]);
			unit = new G4UnionSolid("center+mom+dau2", unit, dBranches[l], rot3, cone_unit_vec3*moving_dist3);
		}
		//construct logical and physical
		if(l==(num_Layers-1)) {
			G4String name = G4UIcommand::ConvertToString(mid_gen+1)+"-gen_"+G4UIcommand::ConvertToString(l)+"-layer";
			G4LogicalVolume* aUnit_log = new G4LogicalVolume(unit, mat_map[146], name+"_log");
            phys_vec.push_back(  new G4PVPlacement(0,                             // Rotation matrix pointer
												   G4ThreeVector(),               // Translation vector
												   aUnit_log ,                    // Logical volume
												   name+"_phy",    			      // Name
												   pre_unit_log,                  // Mother volume
												   false,                         // Unused boolean
												   146+BB1000_bb2000,false));     // Copy number
		}
		else{
            if(168-l+BB1000_bb2000==1167) continue;
            if(168-l+BB1000_bb2000==1166) continue;
            if(168-l+BB1000_bb2000==1161) continue;
            if(168-l+BB1000_bb2000==1160) continue;
            if(168-l+BB1000_bb2000==2167) continue;
            if(168-l+BB1000_bb2000==2166) continue;
            if(168-l+BB1000_bb2000==2163) continue;
            if(168-l+BB1000_bb2000==2162) continue;

			G4String name = G4UIcommand::ConvertToString(mid_gen+1)+"-gen_"+G4UIcommand::ConvertToString(l)+"-layer";
            G4LogicalVolume* aUnit_log = new G4LogicalVolume(unit, mat_map[168-l], name+"_log");
			if(l == 0)
                phys_vec.push_back(  new G4PVPlacement(0,                            // Rotation matrix pointer
													   fBB_vec[1].first - box_center,// Translation vector
													   aUnit_log ,                   // Logical volume
													   name+"_phy",     			 // Name
													   bBoxLogical,                  // Mother volume
													   false,                        // Unused boolean
                                                       168-l+BB1000_bb2000,false));  // Copy number
			else if(l > 0)
                phys_vec.push_back( new G4PVPlacement(0,                             // Rotation matrix pointer
													  G4ThreeVector(),               // Translation vector
													  aUnit_log ,                    // Logical volume
													  name+"_phy",  				 // Name
													  pre_unit_log,                  // Mother volume
													  false,                         // Unused boolean
													  168-l + BB1000_bb2000,false)); // Copy number

			pre_unit_log = aUnit_log;
		}
	}
	return;
}

void LungParallelDetCon::ConstructSD(){
	if(!fSDConstructed){
		fSDConstructed = true;
		SetUpDetectors();
	}
}

void LungParallelDetCon::SetUpDetectors(){
	G4SDManager* pSDman = G4SDManager::GetSDMpointer();
	G4MultiFunctionalDetector* MFDet = new G4MultiFunctionalDetector("lungSD");
	G4VPrimitiveScorer* scorer = new G4PSEnergyDeposit("eDep");

	MFDet->RegisterPrimitive(scorer);

	pSDman->AddNewDetector( MFDet );

    for(auto phyV: phys_vec){
        if(phyV->GetCopyNo()==1146||
           phyV->GetCopyNo()==2146) continue;

        SetSensitiveDetector(phyV->GetLogicalVolume(), MFDet);
	}
}

