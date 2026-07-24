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
// TETModelImport.cc
// \file   MRCP_GEANT4/External/src/TETModelImport.cc
// \author Haegin Han
//

#include "TETModelImport.hh"

namespace {
G4String TrimLocal(G4String str)
{
	const char* whitespace = " \t\r\n";
	size_t begin = str.find_first_not_of(whitespace);
	if(begin == std::string::npos) return "";
	size_t end = str.find_last_not_of(whitespace);
	return str.substr(begin, end - begin + 1);
}

G4String StripNodeExtension(G4String prefix)
{
	const G4String suffix = ".node";
	if(prefix.size() >= suffix.size() &&
	   prefix.substr(prefix.size() - suffix.size()) == suffix){
		return prefix.substr(0, prefix.size() - suffix.size());
	}
	return prefix;
}

G4String StripPhantomDataExtension(G4String prefix)
{
	const std::vector<G4String> suffixes = {".ele", ".material", ".dose", ".RBMnBS", ".DRF", ".hands"};
	for(auto suffix : suffixes){
		if(prefix.size() >= suffix.size() &&
		   prefix.substr(prefix.size() - suffix.size()) == suffix){
			return prefix.substr(0, prefix.size() - suffix.size());
		}
	}
	return prefix;
}

G4String DefaultPhantomDataPrefix(G4String nodePrefix)
{
	const G4String fixedName = "MRCP_AM";
	size_t slash = nodePrefix.find_last_of("/\\");
	if(slash == std::string::npos) return fixedName;
	return nodePrefix.substr(0, slash + 1) + fixedName;
}

void AddUniqueDoseID(std::vector<G4int>& doseIDs, G4int doseID)
{
	if(std::find(doseIDs.begin(), doseIDs.end(), doseID) == doseIDs.end()){
		doseIDs.push_back(doseID);
	}
}
}

TETModelImport::TETModelImport(G4String _phantomName, G4UIExecutive* ui)
:TETModelImport(_phantomName, std::vector<G4String>(), ui)
{}

TETModelImport::TETModelImport(G4String _phantomName, const std::vector<G4String>& extraTetNames, G4UIExecutive* ui)
:TETModelImport(_phantomName, "", extraTetNames, ui)
{}

TETModelImport::TETModelImport(G4String _phantomName, G4String _phantomDataName, const std::vector<G4String>& extraTetNames, G4UIExecutive* ui)
:TETModelImport(_phantomName, _phantomDataName, extraTetNames, std::vector<G4String>(), ui)
{}

TETModelImport::TETModelImport(G4String _phantomName, G4String _phantomDataName, const std::vector<G4String>& extraTetNames, const std::vector<G4String>& patientTetNames, G4UIExecutive* ui)
:doseOrganized(false), hasGeometryBounds(false)
{
	// set phantom name
	phantomName = StripNodeExtension(_phantomName);
	G4String phantomDataPrefix = _phantomDataName.empty()
	                           ? DefaultPhantomDataPrefix(phantomName)
	                           : StripPhantomDataExtension(_phantomDataName);

	G4cout << "================================================================================"<<G4endl;
	G4cout << "\t" << phantomName << " was implemented in this CODE!!   "<< G4endl;
	G4cout << "\tBase phantom node file prefix: " << phantomName << G4endl;
	G4cout << "\tBase phantom data file prefix: " << phantomDataPrefix << G4endl;
	G4cout << "================================================================================"<<G4endl;

	G4String nodeFile     =  phantomName + ".node";
	G4String eleFile      =  phantomDataPrefix + ".ele";
	G4String materialFile =  phantomDataPrefix + ".material";
	G4String doseFile     =  phantomDataPrefix + ".dose";
	G4String boneFile     =  phantomDataPrefix + ".RBMnBS";
	G4String drfFile      =  phantomDataPrefix + ".DRF";
	G4String handsFile    =  phantomDataPrefix + ".hands";

	// read dose file (*.dose) -if there is any
	DoseRead(doseFile);
	// read phantom data files (*.ele, *.node)
    DataRead(eleFile, nodeFile);
	for(auto extraTetName : extraTetNames){
		G4cout << "  Merging additional TET model '" << extraTetName << "'" << G4endl;
		DataRead(extraTetName + ".ele", extraTetName + ".node");
	}
	for(auto patientTetName : patientTetNames){
		G4cout << "  Merging patient TET model '" << patientTetName
		       << "' with material IDs remapped to -(100000 + originalID)" << G4endl;
		DataRead(patientTetName + ".ele", patientTetName + ".node", true);
	}
	FinalizeGeometry();
	// read material file (*.material)
	MaterialRead(materialFile);
	for(auto extraTetName : extraTetNames){
		MaterialRead(extraTetName + ".material");
	}
	for(auto patientTetName : patientTetNames){
		MaterialRead(patientTetName + ".material", true);
	}
	BuildMaterials();
	// read hand element classification file (*.hands) -if there is any
	HandsRead(handsFile);
	// read bone file (*.RBMnBS)
	RBMBSRead(boneFile);
	// read bone file (*.DRF)
	DRFRead(drfFile);
	// read colour data file (colour.dat) if this is interactive mode
	if(ui) ColourRead();
	// print the summary of phantom information
	PrintMaterialInfomation();
}


void TETModelImport::DoseRead(G4String doseFile){
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

void TETModelImport::DataRead(G4String eleFile, G4String nodeFile, G4bool patientNamespace)
{
	G4String tempStr;
	G4int tempInt;
	std::map<G4int, G4int> nodeIndexMap;

	// Read *.node file
	//
	std::ifstream ifpNode;

	ifpNode.open(nodeFile);
	if(!ifpNode.is_open()) {
		// exception for the case when there is no *.node file
		G4Exception("TETModelImport::DataRead","",FatalErrorInArgument,
				G4String("      There is no " + nodeFile ).c_str());
	}
	G4cout << "  Opening TETGEN node (vertex points: x y z) file '" << nodeFile << "'" <<G4endl;

	G4int numVertex;
	G4double xPos, yPos, zPos;

	ifpNode >> numVertex >> tempInt >> tempInt >> tempInt;

	for(G4int i=0; i<numVertex; i++)
	{
		G4int nodeID;
		ifpNode >> nodeID >> xPos >> yPos >> zPos;

		// set the unit
		xPos*=cm;
		yPos*=cm;
		zPos*=cm;

		// save the node data as the form of std::vector<G4ThreeVector>
		nodeIndexMap[nodeID] = vertexVector.size();
		vertexVector.push_back(G4ThreeVector(xPos, yPos, zPos));

		// to get the information of the bounding box of all imported TET models
		if(!hasGeometryBounds){
			boundingBox_Min = G4ThreeVector(xPos, yPos, zPos);
			boundingBox_Max = G4ThreeVector(xPos, yPos, zPos);
			hasGeometryBounds = true;
		}
		else {
			if (xPos < boundingBox_Min.x()) boundingBox_Min.setX(xPos);
			if (xPos > boundingBox_Max.x()) boundingBox_Max.setX(xPos);
			if (yPos < boundingBox_Min.y()) boundingBox_Min.setY(yPos);
			if (yPos > boundingBox_Max.y()) boundingBox_Max.setY(yPos);
			if (zPos < boundingBox_Min.z()) boundingBox_Min.setZ(zPos);
			if (zPos > boundingBox_Max.z()) boundingBox_Max.setZ(zPos);
		}
	}
	ifpNode.close();

	// Read *.ele file
	//
	std::ifstream ifpEle;

	ifpEle.open(eleFile);
	if(!ifpEle.is_open()) {
		// exception for the case when there is no *.ele file
		G4Exception("TETModelImport::DataRead","",FatalErrorInArgument,
				G4String("      There is no " + eleFile ).c_str());
	}
	G4cout << "  Opening TETGEN elements (tetrahedron with node No.) file '" << eleFile << "'" <<G4endl;

	G4int numEle;
	ifpEle >> numEle  >> tempInt >> tempInt;

	for(G4int i=0; i<numEle; i++)
	{
		ifpEle >> tempInt;
		G4int elementID = tempInt;
		G4int* ele = new G4int[4];
		for(G4int j=0;j<4;j++){
			ifpEle >> tempInt;
			if(nodeIndexMap.find(tempInt)==nodeIndexMap.end()){
				G4Exception("TETModelImport::DataRead","",FatalErrorInArgument,
						G4String("      Element in " + eleFile + " references unknown node ID").c_str());
			}
			ele[j]=nodeIndexMap[tempInt];
		}
		G4int copyNo = (G4int)eleVector.size();
		elementIDToCopyNo.insert(std::make_pair(elementID, copyNo));
		elementIDToCopyNos[elementID].push_back(copyNo);
		eleVector.push_back(ele);
		ifpEle >> tempInt;
		materialVector.push_back(patientNamespace ? PatientMaterialID(tempInt) : tempInt);
	}
	ifpEle.close();
}

void TETModelImport::FinalizeGeometry()
{
	if(!hasGeometryBounds) return;

	// Center every imported TET model together in one common bounding box.
	G4ThreeVector center = (boundingBox_Max+boundingBox_Min)*0.5;
	phantomSize = boundingBox_Max - boundingBox_Min;
	std::transform(vertexVector.begin(), vertexVector.end(), vertexVector.begin(),
			[center](G4ThreeVector& v){return v - center;});

	for(size_t i=0; i<eleVector.size(); i++)
	{
		G4int* ele = eleVector[i];
		tetVector.push_back(new G4Tet("Tet_Solid",
									  vertexVector[ele[0]],
									  vertexVector[ele[1]],
									  vertexVector[ele[2]],
									  vertexVector[ele[3]]));

		G4int materialID = materialVector[i];
		if(volumeMap.find(materialID)!=volumeMap.end()){
			volumeMap[materialID] += tetVector[i]->GetCubicVolume();
			numTetMap[materialID]++;
		}
		else {
			volumeMap[materialID] = tetVector[i]->GetCubicVolume();
			numTetMap[materialID] = 1;
		}
	}
}

G4ThreeVector TETModelImport::GetTetrahedronCentroid(G4int idx)
{
	if(idx < 0 || idx >= (G4int)eleVector.size()) return G4ThreeVector();

	G4int* ele = eleVector[idx];
	return (vertexVector[ele[0]]
	      + vertexVector[ele[1]]
	      + vertexVector[ele[2]]
	      + vertexVector[ele[3]]) * 0.25;
}

std::vector<G4int> TETModelImport::GetCopyNumbersForMaterials(const std::set<G4int>& materialIDs)
{
	std::vector<G4int> copyNumbers;
	for(G4int i=0; i<(G4int)materialVector.size(); i++){
		if(materialIDs.find(materialVector[i]) != materialIDs.end()) copyNumbers.push_back(i);
	}
	return copyNumbers;
}

std::vector<G4int> TETModelImport::GetCopyNumbersForElementIDs(const std::set<G4int>& elementIDs,
                                                               const std::set<G4int>& materialIDs)
{
	std::vector<G4int> copyNumbers;
	for(auto elementID : elementIDs){
		auto found = elementIDToCopyNos.find(elementID);
		if(found == elementIDToCopyNos.end()) continue;

		for(auto copyNo : found->second){
			if(copyNo < 0 || copyNo >= (G4int)materialVector.size()) continue;
			if(!materialIDs.empty() && materialIDs.find(materialVector[copyNo]) == materialIDs.end()) continue;
			copyNumbers.push_back(copyNo);
		}
	}
	std::sort(copyNumbers.begin(), copyNumbers.end());
	copyNumbers.erase(std::unique(copyNumbers.begin(), copyNumbers.end()), copyNumbers.end());
	return copyNumbers;
}

G4bool TETModelImport::GetMaterialSetCenterOfMass(const std::set<G4int>& materialIDs, G4ThreeVector& center)
{
	G4double totalVolume = 0.;
	G4ThreeVector weightedCenter;

	for(G4int i=0; i<(G4int)materialVector.size(); i++){
		if(materialIDs.find(materialVector[i]) == materialIDs.end()) continue;
		G4double volume = tetVector[i]->GetCubicVolume();
		weightedCenter += GetTetrahedronCentroid(i) * volume;
		totalVolume += volume;
	}

	if(totalVolume <= 0.) return false;

	center = weightedCenter / totalVolume;
	return true;
}

void TETModelImport::MaterialRead(G4String materialFile, G4bool patientNamespace)
{
	// Read material file (*.material)
	//
	std::ifstream ifpMat;

	ifpMat.open(materialFile);
	if(!ifpMat.is_open()) {
		// exception for the case when there is no *.material file
		G4Exception("TETModelImport::DataRead","",FatalErrorInArgument,
				G4String("      There is no " + materialFile ).c_str());
	}

	G4cout << "  Opening material file '" << materialFile << "'" <<G4endl;

	char read_data[50];
	char* token;
	G4double zaid;
	G4double fraction;
	G4double density;

	while(!ifpMat.eof())
	{
		ifpMat >> read_data;                   //ex) 'C' RBM
		G4String MaterialName;
		ifpMat >> MaterialName;                //ex)  C 'RBM'
		if(MaterialName.empty()) continue;
		ifpMat >> read_data;
		density = std::atof(read_data);        //ex) 1.30
		ifpMat >> read_data;                   //ex) g/cm3
		ifpMat >> read_data;
        if(G4String(read_data).empty()) continue;
		token = std::strtok(read_data,"m");
		G4int matID = std::atoi(token);        //ex) m'10'
		if(patientNamespace) matID = PatientMaterialID(matID);
		G4bool newMaterial = (materialIndexMap.find(matID)==materialIndexMap.end());
		std::vector<std::pair<G4int, G4double>> composition;

		for(G4int i=0 ;  ; i++)
		{
			ifpMat >> read_data;
			if(std::strcmp(read_data, "C")==0 || ifpMat.eof()) break;

			zaid = (G4int)(std::atoi(read_data)/1000);
			ifpMat >> read_data;
			fraction = -1.0 * std::atof(read_data);
			composition.push_back(std::make_pair(G4int(zaid), fraction));
		}

		if(newMaterial){
			materialIndex.push_back(matID);
			organNameMap[matID]= patientNamespace ? "patient_" + MaterialName : MaterialName;
			densityMap[matID] = density*g/cm3;
			materialIndexMap[matID] = composition;
		}
	}
	ifpMat.close();

}

void TETModelImport::BuildMaterials()
{
	// Construct materials for each organ after all material files are read.
	G4Element *elH = new G4Element("TS_H_of_Water", "H", 1., 1.01*g/mole);
	G4NistManager* nistManager = G4NistManager::Instance();

	for(G4int i=0;i<(G4int)materialIndex.size();i++){
		G4int idx = materialIndex[i];
		G4Material* mat = new G4Material(organNameMap[idx], densityMap[idx], G4int(materialIndexMap[idx].size()), kStateSolid, NTP_Temperature, STP_Pressure);
		for(G4int j=0;j<G4int(materialIndexMap[idx].size());j++){
			if(materialIndexMap[idx][j].first==1) mat->AddElement(elH, materialIndexMap[idx][j].second);
			else mat->AddElement(nistManager->FindOrBuildElement(materialIndexMap[idx][j].first), materialIndexMap[idx][j].second);
		}
		materialMap[idx]=mat;
		massMap[idx]=densityMap[idx]*volumeMap[idx];
	}

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

void TETModelImport::HandsRead(G4String handsFile)
{
	std::ifstream ifs(handsFile);
	if(!ifs.is_open()) {
		G4Exception("TETModelImport::HandsRead","",JustWarning,
				G4String("      There is no " + handsFile ).c_str());
		return;
	}

	G4cout << "  Opening hand element classification file '" << handsFile << "'" << G4endl;

	std::set<G4int> leftHand;
	std::set<G4int> rightHand;
	std::set<G4int>* activeSet = nullptr;

	G4String line;
	while(std::getline(ifs, line)){
		size_t comment = line.find('#');
		if(comment != std::string::npos) line = line.substr(0, comment);
		line = TrimLocal(line);
		if(line.empty()) continue;

		if(line == "[left_hand]"){
			activeSet = &leftHand;
			continue;
		}
		if(line == "[right_hand]"){
			activeSet = &rightHand;
			continue;
		}
		if(line.front() == '['){
			activeSet = nullptr;
			continue;
		}
		if(!activeSet) continue;

		std::stringstream ss(line);
		G4String token;
		while(ss >> token){
			size_t dash = token.find('-');
			if(dash == std::string::npos){
				G4int elementID = std::atoi(token.c_str());
				auto copyNo = elementIDToCopyNo.find(elementID);
				if(copyNo != elementIDToCopyNo.end()) activeSet->insert(copyNo->second);
				continue;
			}

			G4int firstID = std::atoi(token.substr(0, dash).c_str());
			G4int lastID = std::atoi(token.substr(dash + 1).c_str());
			if(lastID < firstID) std::swap(firstID, lastID);
			for(G4int elementID = firstID; elementID <= lastID; elementID++){
				auto copyNo = elementIDToCopyNo.find(elementID);
				if(copyNo != elementIDToCopyNo.end()) activeSet->insert(copyNo->second);
			}
		}
	}
	ifs.close();

	handDoseName[LEFT_HAND_DOSE_ID] = "LeftHand";
	handDoseName[RIGHT_HAND_DOSE_ID] = "RightHand";
	handDoseName[BOTH_HANDS_DOSE_ID] = "BothHands";
	handDoseName[LEFT_HAND_SKIN_SENSITIVE_DOSE_ID] = "LeftHand_12401";
	handDoseName[RIGHT_HAND_SKIN_SENSITIVE_DOSE_ID] = "RightHand_12401";
	handDoseName[BOTH_HANDS_SKIN_SENSITIVE_DOSE_ID] = "BothHands_12401";

	for(auto doseNamePair : handDoseName) handDoseMassMap[doseNamePair.first] = 0.;

	auto addHandDose = [this](G4int copyNo, G4int totalDoseID, G4int sensitiveDoseID) {
		if(copyNo < 0 || copyNo >= (G4int)tetVector.size()) return;

		AddUniqueDoseID(handTet2dose[copyNo], totalDoseID);
		AddUniqueDoseID(handTet2dose[copyNo], BOTH_HANDS_DOSE_ID);

		G4int materialID = materialVector[copyNo];
		G4double tetMass = tetVector[copyNo]->GetCubicVolume() * densityMap[materialID];
		handDoseMassMap[totalDoseID] += tetMass;

		if(materialID == HAND_SKIN_SENSITIVE_MATERIAL_ID){
			AddUniqueDoseID(handTet2dose[copyNo], sensitiveDoseID);
			AddUniqueDoseID(handTet2dose[copyNo], BOTH_HANDS_SKIN_SENSITIVE_DOSE_ID);
			handDoseMassMap[sensitiveDoseID] += tetMass;
		}
	};

	for(auto copyNo : leftHand) addHandDose(copyNo, LEFT_HAND_DOSE_ID, LEFT_HAND_SKIN_SENSITIVE_DOSE_ID);
	for(auto copyNo : rightHand) addHandDose(copyNo, RIGHT_HAND_DOSE_ID, RIGHT_HAND_SKIN_SENSITIVE_DOSE_ID);

	std::set<G4int> bothHands(leftHand);
	bothHands.insert(rightHand.begin(), rightHand.end());
	for(auto copyNo : bothHands){
		if(copyNo < 0 || copyNo >= (G4int)tetVector.size()) continue;
		G4int materialID = materialVector[copyNo];
		G4double tetMass = tetVector[copyNo]->GetCubicVolume() * densityMap[materialID];
		handDoseMassMap[BOTH_HANDS_DOSE_ID] += tetMass;
		if(materialID == HAND_SKIN_SENSITIVE_MATERIAL_ID){
			handDoseMassMap[BOTH_HANDS_SKIN_SENSITIVE_DOSE_ID] += tetMass;
		}
	}

	G4cout << "  Hand dose groups: left=" << leftHand.size()
		   << ", right=" << rightHand.size()
		   << ", both=" << bothHands.size() << G4endl;
}

void TETModelImport::RBMBSRead(G4String bonefile){
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
        if(rbmRatio.find(idx)!=rbmRatio.end()) {
            G4cerr<<idx<<" is duplicated in RBMBS file.."<<G4endl;
            exit(0);
        }
		rbmRatio[idx]=rbm;
		bsRatio[idx]=bs;
	}
}

void TETModelImport::DRFRead(G4String DRFfile){
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

void TETModelImport::ColourRead()
{
  // Read colour data file (colour.dat)
  //
  std::ifstream ifpColour;

  ifpColour.open( "colour.dat");
  if(!ifpColour.is_open())
  {
    // exception for the case when there is no colour.dat file
    G4Exception("TETModelImport::DataRead","",FatalErrorInArgument,
                G4String("Colour data file was not found ").c_str());
  }

  G4cout << "  Opening colour data file 'colour.dat'" <<G4endl;

  G4int organID;
  G4double red, green, blue, alpha;
  while( ifpColour >> organID >> red >> green >> blue >> alpha )
  {
    colourMap[organID] = G4Colour(red, green, blue, alpha);
  }

  ifpColour.close();
}

void TETModelImport::PrintMaterialInfomation()
{
	// Print the overall information for each organ
	//
	G4cout << G4endl
		   << std::setw(9)  << "Organ ID"
		   << std::setw(11) << "# of Tet"
		   << std::setw(11) << "vol [cm3]"
		   << std::setw(11) << "d [g/cm3]"
		   << std::setw(11) << "mass [g]"
		   << "\t" << "organ/tissue"<< G4endl ;
	G4cout << "--------------------------------------------------------------------------------"<<G4endl;

	std::map<G4int, G4Material*>::iterator matIter;
	G4cout<<std::setiosflags(std::ios::fixed);
	G4cout.precision(3);
	for(matIter=materialMap.begin(); matIter!=materialMap.end();matIter++)
	{
		G4int idx = matIter->first;

		G4cout << std::setw(9)  << idx                         // organ ID
			   << std::setw(11) << numTetMap[idx]              // # of tetrahedrons
			   << std::setw(11) << volumeMap[idx]/cm3          // organ volume
			   << std::setw(11) << materialMap[idx]
			                       ->GetDensity()/(g/cm3)      // organ density
			   << std::setw(11) << massMap[idx]/g              // organ mass
			   << "\t"<<materialMap[idx]->GetName() << G4endl; // organ name
	}
}
