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

#include "ImportVoxelPhantom.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

ImportVoxelPhantom::ImportVoxelPhantom(G4String _inp)
:inp_name(_inp)
{
	ifstream ifp(inp_name);
	if(!ifp.is_open()) {
		G4Exception("ImportVoxelPhantom::ImportVoxelPhantom","",FatalErrorInArgument,
		G4String("      There is no " + inp_name ).c_str());
	}
	G4String dump1, dump2;

	while (ifp >> dump1)
	{
		if (dump1 == "file=")
		{
			ifp >> dump2;
			lat_name = dump2;
			G4cout << "LATTICE FILE: " << lat_name << G4endl;
			break;
		}
	}
	ifp.close();
		
	ImportPhantomInput();
	ImportPhantomVoxelMaterial();
	ImportPhantomVoxelData();
	ImportPhantomVoxelVolume();
}

ImportVoxelPhantom::~ImportVoxelPhantom() {
	// TODO Auto-generated destructor stub
}


void ImportVoxelPhantom::ImportPhantomInput(){

	ifstream ifp1;
	ifp1.open(lat_name);
	string dump;
	size_t pos;

	// Read Array Size
	getline(ifp1, dump);
	pos = dump.find("0:");
	voxelResolution.push_back(stoi(dump.substr(pos + 2, 4))+1);
	pos = dump.find("0:",pos + 1);
	voxelResolution.push_back(stoi(dump.substr(pos + 2, 4))+1);
	pos = dump.find("0:",pos + 1);
	voxelResolution.push_back(stoi(dump.substr(pos + 2, 4))+1);
		
	G4cout << "VOX RES: " << voxelResolution[0] << " " <<voxelResolution[1] << " " <<voxelResolution[2] << G4endl;
	
	ifp1.close();
	
	// Read VoxelSize
	ifstream ifp2;
	ifp2.open(inp_name);
	string dump1, dump2;
	G4cout<<inp_name<<G4endl;
	while (!ifp2.eof()) {
	ifp2 >> dump1;
	if (dump1 == "110") {
		ifp2 >> dump2;
		if (dump2 == "rpp") {
			ifp2 >> dump1;
			ifp2 >> dump1;
			voxelSize.setX(stod(dump1)*cm);
			ifp2 >> dump1;
			ifp2 >> dump1;
			voxelSize.setY(stod(dump1)*cm);
			ifp2 >> dump1;
			ifp2 >> dump1;
			voxelSize.setZ(stod(dump1)*cm);
			break;
		}
	}
	
	}
	ifp2.close();
	//cout << voxelSize << endl;

	phantomSize.setX(voxelResolution[0]*voxelSize.x());
	phantomSize.setY(voxelResolution[1]*voxelSize.y());
	phantomSize.setZ(voxelResolution[2]*voxelSize.z());
	G4cout << "VOX SIZE[mm]: " << phantomSize << G4endl;

	// Read DAP Volume
	ifstream ifp3;
	ifp3.open(inp_name);

	while (!ifp3.eof())
	{
	ifp3 >> dump1;
	if (dump1 == "1062")
	{
		ifp3 >> dump2;
		ifp3 >> dump2;
		ifp3 >> dump2;
		DAPHalfSizeX = stod(dump2) * cm;
		G4cout << "DAPHalfSizeX: " << DAPHalfSizeX << G4endl;
	}

	else if (dump1 == "2062")
	{
		ifp3 >> dump2;
		ifp3 >> dump2;
		ifp3 >> dump2;
		DAPHalfSizeY = stod(dump2) * cm;
		G4cout << "DAPHalfSizeY: " << DAPHalfSizeY << G4endl;
	}

	else if (dump1 == "3062")
	{
		ifp3 >> dump2;
		ifp3 >> dump2;
		ifp3 >> dump2;
		DAPHalfSizeZ = stod(dump2) * cm;
		G4cout << "DAPHalfSizeZ: " << DAPHalfSizeZ << G4endl;
	}
	}
	ifp3.close();	
}

void ImportVoxelPhantom::ImportPhantomVoxelData(){

		ifstream ifp1;
		ifp1.open(lat_name);;
		string dump;
		vector<string> strArr;
		size_t pos;

		// Read Phantom Voxel Array Data
		while (!ifp1.eof()) {
		//for (int k=0; k<2; k++){
			getline(ifp1, dump);	
			pos = dump.find(" ");
				vector<size_t> posArr;
				posArr.push_back(0);
				posArr.push_back(pos);
				while (pos != string::npos) {
					//cout << pos << endl;
					pos = dump.find(" ", pos + 1);
					posArr.push_back(pos);
				}
				
				if (posArr.size() > 5) {
					size_t pos_i = posArr[0] + posArr[1] + posArr[2] + posArr[3] + posArr[4];
					//cout << posArr[0] << " " << posArr[1] << " " << posArr[2] << " " << posArr[3] << " " << posArr[4] << " " << endl;
					//size_t pos_i = 1;
					//cout << pos_i << endl;
					if (pos_i == 6) {
						//strArr.push_back(dump.substr(posArr[0], posArr[1] - posArr[0]));
						//cout << dump.substr(posArr[0], posArr[1] - posArr[0]) << endl;
						for (int i = 2; i < posArr.size(); i++)
						{
							if ((posArr[i] - posArr[i - 1]) > 1) {
								strArr.push_back(dump.substr(posArr[i - 1] + 1, posArr[i] - posArr[i - 1] - 1));
								//cout << dump.substr(posArr[i - 1] + 1, posArr[i] - posArr[i - 1] - 1) << endl;
							}
						}
						
					}
				}
				posArr.clear();
		}
		ifp1.close();
		
	//cout << strArr[0] << endl;
	//cout << strArr[strArr.size()-1] << endl;
	//pos = strArr[strArr.size() - 1].find("r");
	//cout << stoi(strArr[strArr.size() - 1].substr(0, pos)) << endl;
	
	voxelData = new G4int**[voxelResolution[0]];
	for(int i=0;i<voxelResolution[0];i++){
		voxelData[i] = new G4int*[voxelResolution[1]];
		for(int j=0;j<voxelResolution[1];j++){
			voxelData[i][j] = new G4int[voxelResolution[2]];
			for(int k=0;k<voxelResolution[2];k++){
				voxelData[i][j][k] = 0;
			}
		}
	}
		
	int VoxelCount = 0;
	int VoxelDump = 0;
	for (int kk = 0; kk < strArr.size(); kk++)
	{
		pos = strArr[kk].find("r");
		if (pos != string::npos) {			
			int r = stoi(strArr[kk].substr(0, pos));
			for (int jj = 0; jj < r; jj++) {
				VoxelDump = stoi(strArr[kk - 1]);
				int x = (VoxelCount % (voxelResolution[0] * voxelResolution[1])) % voxelResolution[0];
				int y = (VoxelCount % (voxelResolution[0] * voxelResolution[1])) / voxelResolution[0];
				int z = VoxelCount / (voxelResolution[0] * voxelResolution[1]);
				voxelData[x][y][z] = VoxelDump;
				numVoxel[voxelData[x][y][z]]++;
				VoxelCount = VoxelCount + 1;
			}
		}
		else {
			VoxelDump = stoi(strArr[kk]);
			int x = (VoxelCount % (voxelResolution[0] * voxelResolution[1])) % voxelResolution[0];
			int y = (VoxelCount % (voxelResolution[0] * voxelResolution[1])) / voxelResolution[0];
			int z = VoxelCount / (voxelResolution[0] * voxelResolution[1]);
			voxelData[x][y][z] = VoxelDump;
			numVoxel[voxelData[x][y][z]]++;
			VoxelCount = VoxelCount + 1;
			
		}
	}
	strArr.clear();

		int TotalNumVoxel(0);
		for(int k=0; k<255; k++)
		{
			G4cout<<k<<":   "<<numVoxel[k]<<G4endl;
		}
			G4cout<<"totalNumVoxel_check1:   "<<VoxelCount<<G4endl;
			G4cout<<"totalNumVoxel_check2:   "<<voxelResolution[0]*voxelResolution[1]*voxelResolution[2]<<G4endl;
}

void ImportVoxelPhantom::ImportPhantomVoxelMaterial(){


	for (int i=0; i<256; i++)
	{
		CompID[i] = 0; //VoxID -> CompID
		Density[i] = 0; //VoxID -> Density
	}
	
	size_t pos;
	
	ifstream ifp3;
	ifp3.open(lat_name);
	string dump;
	getline(ifp3, dump);
	while (!ifp3.eof()) {
		//for (int k=0; k<2; k++){
		getline(ifp3, dump);
		pos = dump.find(" ");
		vector<size_t> posArr;
		posArr.push_back(0);
		posArr.push_back(pos);
		while (pos != string::npos) {
			//cout << pos << endl;
			pos = dump.find(" ", pos + 1);
			posArr.push_back(pos);
		}

		if (posArr.size() > 5) {
			size_t pos_i = posArr[0] + posArr[1] + posArr[2] + posArr[3] + posArr[4];
			//cout << posArr[0] << " " << posArr[1] << " " << posArr[2] << " " << posArr[3] << " " << posArr[4] << " " << endl;
			//size_t pos_i = 1;
			//cout << pos_i << endl;
			if (pos_i != 6) {
				CompID[stoi(dump.substr(0, 3))] = stoi(dump.substr(9, 2));
				Density[stoi(dump.substr(0, 3))] = stod(dump.substr(14, 5));
				//cout << dump.substr(0, 3) << " " << CompID[stoi(dump.substr(0, 3))] << " " << Density[stoi(dump.substr(0, 3))] << endl;
			}
		}
		posArr.clear();
	}
	ifp3.close();

	// Read Composition
	ifstream ifp4;
	ifp4.open(inp_name);
	map<int, vector<int>> Zaid;
	map<int, vector<double>> Fraction;
	
	int CompID_dump = 0;
	
	int Cn = 0;
	while (!ifp4.eof()) {
		//getline(ifp4, dump);
		ifp4 >> dump;
		//cout << dump.size() << endl;
		if (dump.size()>0 && dump.at(0) == 'm') {
			if (dump.substr(1, 3).find_first_not_of("0123456789 ") == string::npos) {
				CompID_dump = stoi(dump.substr(1, 3));
				Cn = 0;
				//cout << dump.substr(1, 3) << endl;
			}
		}

		if (CompID_dump != 0) {
			if (dump.find_first_not_of("0123456789-.") == string::npos) {
				Cn++;
				if (Cn % 2 == 1) Zaid[CompID_dump].push_back(stoi(dump)*0.001) ;
				if (Cn % 2 == 0) Fraction[CompID_dump].push_back(abs(stod(dump))*0.01);
			}
		}

		if (dump == "mode") break;

	}
	ifp4.close();
	
	int z = 10;
	for (int i = 0; i < Zaid[z].size(); i++)
	{
		//cout << Zaid[z][i] <<" "<< Fraction[z][i] << endl;
	}

    G4Element *elH = new G4Element("TS_H_of_Water", "H", 1., 1.01*g/mole);
    G4NistManager* nistManager = G4NistManager::Instance();
//
    
	for(G4int VoxID=0;VoxID<256;VoxID++){
		if (CompID[VoxID]!=0){
			
		G4Material* mat = new G4Material(to_string(VoxID), Density[VoxID] * g/cm3, G4int(Zaid[CompID[VoxID]].size()), kStateSolid, NTP_Temperature, STP_Pressure);
		for(G4int j=0;j<G4int(Zaid[CompID[VoxID]].size());j++){
			if(Zaid[CompID[VoxID]][j]==1) mat->AddElement(elH, Fraction[CompID[VoxID]][j]);
            else mat->AddElement(nistManager->FindOrBuildElement(Zaid[CompID[VoxID]][j]), Fraction[CompID[VoxID]][j]);
		}
	    //G4cout<<VoxID<<": "<<mat->GetName()<<" density: "<<mat->GetDensity()/(g/cm3)<<G4endl;
		materialMap[VoxID]=mat;
		materialIndex.push_back(VoxID);
		}
	}
	cout<<"Number of Materials: "<< materialMap.size()<<endl;
}

void ImportVoxelPhantom::ImportPhantomVoxelVolume(){
	for(int k=0; k<256; k++)
	{
		organVolume[k] = numVoxel[k] *  voxelSize[0] * voxelSize[1] * voxelSize[2];
	}
}

