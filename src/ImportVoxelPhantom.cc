/*
 * ImportVoxelPhantom.cc
 *
 *  Created on: May 20, 2016
 *      Author: mchan
 */

#include "ImportVoxelPhantom.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

ImportVoxelPhantom::ImportVoxelPhantom(G4String phantomFile)
:Filename(phantomFile) {

	G4String InformFile = "Phantom_information";
	G4String voxelFile = "./phantoms/"+Filename+"_data";
	G4String materialFile = "./phantoms/"+Filename+"_material";
	G4String DRFfile = "./phantoms/"+Filename+"_DRF";
	G4String RBMBSfile = "./phantoms/"+Filename+"_RBMnBS";


	ImportPhantomInput(InformFile);
	ImportPhantomVoxelMaterial(materialFile);
	ImportPhantomVoxelData(voxelFile);
	ImportPhantomVoxelVolume();
}

ImportVoxelPhantom::~ImportVoxelPhantom() {
}


void ImportVoxelPhantom::ImportPhantomInput(G4String inputFile){

	std::ifstream ifp;
	ifp.open(inputFile.c_str());

	if(!ifp.is_open()) {
		G4cerr << inputFile << " not found!!" << G4endl;
		return;
	}

	char read_data[50];
	G4int fNx, fNy, fNz;
	  while(!ifp.eof())
	  {
		  ifp >> read_data;
		  if( read_data==Filename ) {
			  ifp >> read_data;
			  ifp >> fNx >> fNy >> fNz;
			  voxelResolution.push_back(fNx);
			  voxelResolution.push_back(fNy);
			  voxelResolution.push_back(fNz);
			  ifp >> read_data;
			  ifp >> read_data;
			  voxelSize.setX(std::atof(read_data));
			  ifp >> read_data;
			  voxelSize.setY(std::atof(read_data));
			  ifp >> read_data;
			  voxelSize.setZ(std::atof(read_data));

			  phantomSize.setX(fNx*voxelSize.x());
			  phantomSize.setY(fNy*voxelSize.y());
			  phantomSize.setZ(fNz*voxelSize.z());
			  break;
		  }
	  }

	ifp.close();
}

void ImportVoxelPhantom::ImportPhantomVoxelData(G4String voxelFile){

	std::ifstream ifp;
	ifp.open(voxelFile.c_str());

	if(!ifp.is_open()) {
		G4cerr << voxelFile << " not found!!" << G4endl;
		return;
	}

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


	for(int k=0;k<voxelResolution[2];k++){
		for(int j=0;j<voxelResolution[1];j++){
			for(int i=0;i<voxelResolution[0];i++){
				ifp >> voxelData[i][j][k];
				numVoxel[voxelData[i][j][k]]++;
			}
		}
	}

	ifp.close();

}

void ImportVoxelPhantom::ImportPhantomVoxelMaterial(G4String materialFile){

	std::ifstream ifp;
	ifp.open(materialFile.c_str());

	if(!ifp.is_open()) {
		G4cerr << materialFile << " not found!!" << G4endl;
		return;
	}

	char read_data[50];
	char* token;
	G4double zaid;
	G4double fraction;
	G4String MaterialName;
	G4double density;
	 while(!ifp.eof())
	  {
		 ifp >> read_data;                   //ex) 'C' RBM
		 ifp >> MaterialName;                //ex)  C 'RBM'
		 ifp >> read_data;
		 density = std::atof(read_data);  //ex) 1.30
		 ifp >> read_data;                   //ex) g/cm3
		 ifp >> read_data;
		token = std::strtok(read_data,"m");
		G4int matID = std::atoi(token);                 //ex) m'10'
		materialIndex.push_back(matID);
		organNameMap[matID]= MaterialName;
		densityMap[matID] = density;

		for(G4int i=0 ;  ; i++)
		{
			ifp >> read_data;
			if(std::strcmp(read_data, "C")==0 || ifp.eof()) break;

			zaid = (G4int)(std::atoi(read_data)/1000);
			ifp >> read_data;
			fraction = -1.0 * std::atof(read_data);
			materialIndexMap[matID].push_back(std::make_pair(G4int(zaid), fraction));
		}

	  }

    G4Element *elH = new G4Element("TS_H_of_Water", "H", 1., 1.01*g/mole);
    G4NistManager* nistManager = G4NistManager::Instance();

	for(G4int i=0;i<(G4int)materialIndex.size();i++){
		G4int idx = materialIndex[i];
		G4Material* mat = new G4Material(organNameMap[idx], densityMap[idx] * g/cm3, G4int(materialIndexMap[idx].size()), kStateSolid, NTP_Temperature, STP_Pressure);
		for(G4int j=0;j<G4int(materialIndexMap[idx].size());j++){
			if(materialIndexMap[idx][j].first==1) mat->AddElement(elH, materialIndexMap[idx][j].second);
            else mat->AddElement(nistManager->FindOrBuildElement(materialIndexMap[idx][j].first), materialIndexMap[idx][j].second);
		}
		G4cout<<idx<<": "<<mat->GetName()<<" density: "<<mat->GetDensity()/(g/cm3)<<G4endl;
		materialMap[idx]=mat;
	}

	ifp.close();
}



void ImportVoxelPhantom::ImportPhantomVoxelVolume(){
	for(auto iter:materialIndex)
		organVolume[iter] = numVoxel[iter] *  voxelSize[0] * voxelSize[1] * voxelSize[2];
}

