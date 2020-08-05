#include "VOXModelImport.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

VOXModelImport::VOXModelImport(G4String phantomFile)
:Filename(phantomFile) {

	G4String InformFile = "Phantom_information";
	G4String voxelFile = Filename+"_data";
	G4String materialFile = Filename+"_material";
	G4String doseFile = Filename + ".dose";
	G4String DRFfile = Filename+"_DRF";
	G4String RBMBSfile = Filename+"_RBMnBS";


	ImportPhantomInfo(InformFile);
	ImportPhantomVoxelMaterial(materialFile);
	ImportPhantomVoxelData(voxelFile);
	DoseRead(doseFile);
	ImportPhantomVoxelVolume();
	DRFRead(DRFfile);
	RBMBSRead(RBMBSfile);
}

VOXModelImport::~VOXModelImport() {
}


void VOXModelImport::ImportPhantomInfo(G4String inputFile){

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
		  if( Filename.find(G4String(read_data))!= std::string::npos ) {
			  ifp >> read_data;
			  ifp >> fNx >> fNy >> fNz;
			  voxelResolution.push_back(fNx);
			  voxelResolution.push_back(fNy);
			  voxelResolution.push_back(fNz);
			  ifp >> read_data;
			  ifp >> read_data;
			  voxelSize.setX(std::atof(read_data)*mm);
			  ifp >> read_data;
			  voxelSize.setY(std::atof(read_data)*mm);
			  ifp >> read_data;
			  voxelSize.setZ(std::atof(read_data)*mm);

			  phantomSize.setX(fNx*voxelSize.x());
			  phantomSize.setY(fNy*voxelSize.y());
			  phantomSize.setZ(fNz*voxelSize.z());
			  break;
		  }
	  }

	ifp.close();
}

void VOXModelImport::ImportPhantomVoxelData(G4String voxelFile){

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

void VOXModelImport::ImportPhantomVoxelMaterial(G4String materialFile){

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
		organID.push_back(matID);
		organNameMap[matID]= MaterialName;
		densityMap[matID] = density* g/cm3;

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

	for(G4int i=0;i<(G4int)organID.size();i++){
		G4int idx = organID[i];
		G4Material* mat = new G4Material(organNameMap[idx], densityMap[idx], G4int(materialIndexMap[idx].size()), kStateSolid, NTP_Temperature, STP_Pressure);
		for(G4int j=0;j<G4int(materialIndexMap[idx].size());j++){
			if(materialIndexMap[idx][j].first==1) mat->AddElement(elH, materialIndexMap[idx][j].second);
            else mat->AddElement(nistManager->FindOrBuildElement(materialIndexMap[idx][j].first), materialIndexMap[idx][j].second);
		}
		G4cout<<idx<<": "<<mat->GetName()<<" density: "<<mat->GetDensity()/(g/cm3)<<G4endl;
		materialMap[idx]=mat;
	}
	ifp.close();
}

void VOXModelImport::ImportPhantomVoxelVolume(){
	for(auto iter:organID)
		organVolume[iter] = numVoxel[iter] *  voxelSize[0] * voxelSize[1] * voxelSize[2];

	for(auto mat:organVolume)
		massMap[mat.first] = mat.second * densityMap[mat.first];

	if(DoseWasOrganized()){
		for(auto dm:doseName){
			doseMassMap[dm.first] = 0;
		}
		for(auto od:organ2dose){
			for(auto doseID:od.second)
				doseMassMap[doseID] += massMap[od.first];
		}
	}

}

void VOXModelImport::DoseRead(G4String doseFile){
	//read dose file : PLEASE be careful not to include dose ID 0
	std::ifstream ifs(doseFile);
	if(!ifs.is_open()) return;
	doseOrganized = true;

	G4String aLine;
	while(!ifs.eof()){
		getline(ifs, aLine);
		if(aLine.empty()) break;

		std::stringstream ss(aLine);
		G4int doseID; ss>>doseID;
		G4String name; ss>>name; doseName[doseID] = name;
		G4int organID;
		while(ss>>organID){
			if(organ2dose.find(organID)==organ2dose.end()) organ2dose[organID] = {doseID};
			else	                                       organ2dose[organID].push_back(doseID);
		}
	}
	ifs.close();

}

void VOXModelImport::DRFRead(G4String DRFfile){
	std::ifstream ifp;
	ifp.open(DRFfile.c_str());

	if(!ifp.is_open()) {
		G4cerr << DRFfile << " not found!!" << G4endl;
		return;
	}

	G4int ID;
    G4double DRF;

    while (!ifp.eof()) {
        G4String dump;
        getline(ifp, dump);
        std::stringstream ss(dump); dump.clear();
        ss >> dump;
        if(dump.empty()) continue;
        ID = atoi(dump.c_str());
        if(rbmDRF.find(ID)!=rbmDRF.end()) {
            G4cerr<<ID<<" is duplicated in DRF file.."<<G4endl;
            exit(0);
        }
        rbmDRF[ID]={};
        bsDRF[ID]={};
        for (int j=0; j<25; j++) {
            ss >> DRF;
            rbmDRF[ID].push_back(DRF);
        }
        for (int j=0; j<25; j++) {
            ss >> DRF;
            bsDRF[ID].push_back(DRF);
        }
    }

    ifp.close();
}

void VOXModelImport::RBMBSRead(G4String bonefile){
	std::ifstream ifs(bonefile);
	if(!ifs.is_open()) {
		// exception for the case when there is no *.material file
		G4Exception("TETModelImport::RBMBSRead","",JustWarning,
				G4String("      There is no " + bonefile ).c_str());
		return;
	}
	G4int idx;
	G4double rbm, bs;
	while(ifs>>idx>>rbm>>bs){
		rbmRatio[idx]=rbm;
		bsRatio[idx]=bs;
	}
}

