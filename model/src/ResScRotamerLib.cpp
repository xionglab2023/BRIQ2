/*
 * ResScRotamerLib.cpp
 *
 */

#include "model/ResScRotamerLib.h"

namespace NSPmodel {

ResScRotamerLib::ResScRotamerLib(){
	// TODO Auto-generated constructor stub

	this->atLib = new AtomLib();
	string path = NSPdataio::datapath();
	ResName rn;
	ifstream file;
	string s;
	string fileName;
	vector<string> spt;
	int index;
	char xx[200];

	for(int i=0;i<20;i++){
		string tri = rn.intToTri(i);
		if(tri == "GLY") {
			this->rotList[5].push_back(new ResScRotamer());
			continue;
		}
		fileName = path + "scRotamer/rot/rotLib/" + tri +".rot";
		file.open(fileName.c_str(), ios::in);
		if(! file.is_open()) {
			cout << "fail to open sidechain rotamer library file " << fileName << endl;
			exit(1);
		}

		getline(file, s);
		int index = 0;
		while(getline(file, s)){
			ResScRotamer* rot = new ResScRotamer(s, atLib);
			rot->rotID = index;
			this->rotList[i].push_back(rot);
			index++;
		}
		file.close();
	}

	for(int i=0;i<20;i++){
		this->rotNum.push_back(this->rotList[i].size());
	}

	int idA, idB;
	double ene;
	string tag = "";

	fileName = path + "scRotamer/ene" + tag + "/ALA.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eALA[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/CYS.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eCYS[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/ASP.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eASP[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/GLU.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eGLU[idA][idB] = ene;
	}
	file.close();


	fileName = path + "scRotamer/ene" + tag + "/PHE.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->ePHE[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/HIS.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eHIS[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/ILE.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eILE[idA][idB] = ene;
	}
	file.close();


	fileName = path + "scRotamer/ene" + tag + "/LYS.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eLYS[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/LEU.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eLEU[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/MET.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eMET[idA][idB] = ene;
	}
	file.close();


	fileName = path + "scRotamer/ene" + tag + "/ASN.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eASN[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/PRO.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->ePRO[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/GLN.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eGLN[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/ARG.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eARG[idA][idB] = ene;
	}
	file.close();


	fileName = path + "scRotamer/ene" + tag + "/SER.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eSER[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/THR.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eTHR[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/VAL.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eVAL[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/TRP.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eTRP[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/TYR.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eTYR[idA][idB] = ene;
	}
	file.close();

	/*
	 * GLU, GLN, ARG, LYS, MET cluster: coverage*1
	 */
	fileName = path + "scRotamer/rot/twoLevel/GLU.cluster-x1";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer cluster1K file " << fileName << endl;
		exit(1);
	}
	index = 0;
	while(getline(file, s)){
		splitString(s, " ", &spt);
		int centerRotID = atoi(spt[0].c_str());
		for(int i=0;i<spt.size();i++){
			int rotID = atoi(spt[i].c_str());
			GLUIndex[rotID] = index;
			GLUCenter[rotID] = centerRotID;
		}
		index ++;
	}
	file.close();

	fileName = path + "scRotamer/rot/twoLevel/GLN.cluster-x1";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer cluster1K file " << fileName << endl;
		exit(1);
	}
	index = 0;
	while(getline(file, s)){
		splitString(s, " ", &spt);
		int centerRotID = atoi(spt[0].c_str());
		for(int i=0;i<spt.size();i++){
			int rotID = atoi(spt[i].c_str());
			GLNIndex[rotID] = index;
			GLNCenter[rotID] = centerRotID;
		}
		index ++;
	}
	file.close();

	fileName = path + "scRotamer/rot/twoLevel/LYS.cluster-x1";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer cluster1K file " << fileName << endl;
		exit(1);
	}
	index = 0;
	while(getline(file, s)){
		splitString(s, " ", &spt);
		int centerRotID = atoi(spt[0].c_str());
		for(int i=0;i<spt.size();i++){
			int rotID = atoi(spt[i].c_str());
			LYSIndex[rotID] = index;
			LYSCenter[rotID] = centerRotID;
		}
		index ++;
	}
	file.close();

	fileName = path + "scRotamer/rot/twoLevel/MET.cluster-x1";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer cluster1K file " << fileName << endl;
		exit(1);
	}
	index = 0;
	while(getline(file, s)){
		splitString(s, " ", &spt);
		int centerRotID = atoi(spt[0].c_str());
		for(int i=0;i<spt.size();i++){
			int rotID = atoi(spt[i].c_str());
			METIndex[rotID] = index;
			METCenter[rotID] = centerRotID;
		}
		index ++;
	}
	file.close();

	fileName = path + "scRotamer/rot/twoLevel/ARG.cluster-x1";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer cluster1K file " << fileName << endl;
		exit(1);
	}
	index = 0;
	while(getline(file, s)){
		splitString(s, " ", &spt);
		int centerRotID = atoi(spt[0].c_str());
		for(int i=0;i<spt.size();i++){
			int rotID = atoi(spt[i].c_str());
			ARGIndex[rotID] = index;
			ARGCenter[rotID] = centerRotID;
		}
		index ++;
	}
	file.close();

	/*
	 * GLU, GLN, ARG, LYS, MET local energy change
	 */
	fileName = path + "scRotamer/rot/den/GLU.den";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer local density file " << fileName << endl;
		exit(1);
	}
	double GLUDensity[10000];
	index = 0;
	while(getline(file, s)){
		GLUDensity[index] = log(atof(s.c_str()));
		index++;
	}
	for(int i=0;i<10000;i++){
		GLUELocal[i] = GLUDensity[GLUCenter[i]] - GLUDensity[i];
	}
	file.close();

	fileName = path + "scRotamer/rot/den/GLN.den";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer local density file " << fileName << endl;
		exit(1);
	}
	double GLNDensity[10000];
	index = 0;
	while(getline(file, s)){
		GLNDensity[index] = log(atof(s.c_str()));
		index++;
	}
	for(int i=0;i<10000;i++){
		GLNELocal[i] = GLNDensity[GLNCenter[i]] - GLNDensity[i];
	}
	file.close();

	fileName = path + "scRotamer/rot/den/LYS.den";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer local density file " << fileName << endl;
		exit(1);
	}
	double LYSDensity[10000];
	index = 0;
	while(getline(file, s)){
		LYSDensity[index] = log(atof(s.c_str()));
		index++;
	}
	for(int i=0;i<10000;i++){
		LYSELocal[i] = LYSDensity[LYSCenter[i]] - LYSDensity[i];
	}
	file.close();

	fileName = path + "scRotamer/rot/den/MET.den";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer local density file " << fileName << endl;
		exit(1);
	}
	double METDensity[10000];
	index = 0;
	while(getline(file, s)){
		METDensity[index] = log(atof(s.c_str()));
		index++;
	}
	for(int i=0;i<10000;i++){
		METELocal[i] = METDensity[METCenter[i]] - METDensity[i];
	}
	file.close();

	fileName = path + "scRotamer/rot/den/ARG.den";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer local density file " << fileName << endl;
		exit(1);
	}
	double ARGDensity[20000];
	index = 0;
	while(getline(file, s)){
		ARGDensity[index] = log(atof(s.c_str()));
		index++;
	}
	for(int i=0;i<20000;i++){
		ARGELocal[i] = ARGDensity[ARGCenter[i]] - ARGDensity[i];
	}
	file.close();


	/*
	 * GLU, GLN, ARG, LYS, MET cluster: coverage*4
	 */

	for(int i=0;i<20;i++){
		for(int j=0;j<1000;j++){
			for(int k=0;k<210;k++){
				rotCluster[i][j][k] = 0;
			}
		}
	}

	vector<string> triList;
	triList.push_back("GLU");
	triList.push_back("GLN");
	triList.push_back("LYS");
	triList.push_back("MET");
	triList.push_back("ARG");
	for(int i=0;i<triList.size();i++){
		string tri = triList[i];
		fileName = path + "scRotamer/rot/twoLevel/"+tri+".cluster-x4";
		file.open(fileName.c_str(), ios::in);
		if(! file.is_open()) {
			cout << "fail to open sidechain rotamer cluster file " << fileName << endl;
			exit(1);
		}
		index = 0;
		while(getline(file, s)){
			splitString(s, " ", &spt);
			int centerRotID = atoi(spt[0].c_str());
			rotCluster[rn.triToInt(tri)][index][0] = spt.size();
			for(int i=0;i<spt.size();i++){
				int rotID = atoi(spt[i].c_str());
				rotCluster[rn.triToInt(tri)][index][i+1] = rotID;
			}
			index ++;
		}
		file.close();

        fileName = path + "scRotamer/rot/twoLevel/"+tri+".cluster-x1";
        file.open(fileName.c_str(), ios::in);
        if(! file.is_open()) {
                cout << "fail to open sidechain rotamer cluster file " << fileName << endl;
                exit(1);
        }
        index = 0;
        while(getline(file, s)){
                splitString(s, " ", &spt);
                int centerRotID = atoi(spt[0].c_str());
                rotClusterUnique[rn.triToInt(tri)][index][0] = spt.size();
                for(int i=0;i<spt.size();i++){
                        int rotID = atoi(spt[i].c_str());
                        rotClusterUnique[rn.triToInt(tri)][index][i+1] = rotID;
                }
                index ++;
        }
        file.close();
	}
}

ResScRotamerLib::ResScRotamerLib(const string& tag){
	// TODO Auto-generated constructor stub

	this->atLib = new AtomLib();
	string path = NSPdataio::datapath();
	ResName rn;
	ifstream file;
	string s;
	string fileName;
	vector<string> spt;
	int index;
	char xx[200];

	for(int i=0;i<20;i++){
		string tri = rn.intToTri(i);
		if(tri == "GLY") {
			this->rotList[5].push_back(new ResScRotamer());
			continue;
		}
		fileName = path + "scRotamer/rot/rotLib/" + tri +".rot";
		file.open(fileName.c_str(), ios::in);
		if(! file.is_open()) {
			cout << "fail to open sidechain rotamer library file " << fileName << endl;
			exit(1);
		}

		getline(file, s);
		int index = 0;
		while(getline(file, s)){
			ResScRotamer* rot = new ResScRotamer(s, atLib);
			rot->rotID = index;
			this->rotList[i].push_back(rot);
			index++;
		}
		file.close();
	}

	for(int i=0;i<20;i++){
		this->rotNum.push_back(this->rotList[i].size());
	}

	int idA, idB;
	double ene;

	fileName = path + "scRotamer/ene" + tag + "/ALA.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eALA[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/CYS.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eCYS[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/ASP.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eASP[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/GLU.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eGLU[idA][idB] = ene;
	}
	file.close();


	fileName = path + "scRotamer/ene" + tag + "/PHE.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->ePHE[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/HIS.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eHIS[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/ILE.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eILE[idA][idB] = ene;
	}
	file.close();


	fileName = path + "scRotamer/ene" + tag + "/LYS.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eLYS[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/LEU.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eLEU[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/MET.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eMET[idA][idB] = ene;
	}
	file.close();


	fileName = path + "scRotamer/ene" + tag + "/ASN.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eASN[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/PRO.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->ePRO[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/GLN.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eGLN[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/ARG.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eARG[idA][idB] = ene;
	}
	file.close();


	fileName = path + "scRotamer/ene" + tag + "/SER.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eSER[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/THR.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eTHR[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/VAL.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eVAL[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/TRP.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eTRP[idA][idB] = ene;
	}
	file.close();

	fileName = path + "scRotamer/ene" + tag + "/TYR.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer energy file " << fileName << endl;
		exit(1);
	}
	while(file >> idA >> idB >> ene){
		this->eTYR[idA][idB] = ene;
	}
	file.close();

	/*
	 * GLU, GLN, ARG, LYS, MET cluster: coverage*1
	 */
	fileName = path + "scRotamer/rot/twoLevel/GLU.cluster-x1";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer cluster1K file " << fileName << endl;
		exit(1);
	}
	index = 0;
	while(getline(file, s)){
		splitString(s, " ", &spt);
		int centerRotID = atoi(spt[0].c_str());
		for(int i=0;i<spt.size();i++){
			int rotID = atoi(spt[i].c_str());
			GLUIndex[rotID] = index;
			GLUCenter[rotID] = centerRotID;
		}
		index ++;
	}
	file.close();

	fileName = path + "scRotamer/rot/twoLevel/GLN.cluster-x1";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer cluster1K file " << fileName << endl;
		exit(1);
	}
	index = 0;
	while(getline(file, s)){
		splitString(s, " ", &spt);
		int centerRotID = atoi(spt[0].c_str());
		for(int i=0;i<spt.size();i++){
			int rotID = atoi(spt[i].c_str());
			GLNIndex[rotID] = index;
			GLNCenter[rotID] = centerRotID;
		}
		index ++;
	}
	file.close();

	fileName = path + "scRotamer/rot/twoLevel/LYS.cluster-x1";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer cluster1K file " << fileName << endl;
		exit(1);
	}
	index = 0;
	while(getline(file, s)){
		splitString(s, " ", &spt);
		int centerRotID = atoi(spt[0].c_str());
		for(int i=0;i<spt.size();i++){
			int rotID = atoi(spt[i].c_str());
			LYSIndex[rotID] = index;
			LYSCenter[rotID] = centerRotID;
		}
		index ++;
	}
	file.close();

	fileName = path + "scRotamer/rot/twoLevel/MET.cluster-x1";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer cluster1K file " << fileName << endl;
		exit(1);
	}
	index = 0;
	while(getline(file, s)){
		splitString(s, " ", &spt);
		int centerRotID = atoi(spt[0].c_str());
		for(int i=0;i<spt.size();i++){
			int rotID = atoi(spt[i].c_str());
			METIndex[rotID] = index;
			METCenter[rotID] = centerRotID;
		}
		index ++;
	}
	file.close();

	fileName = path + "scRotamer/rot/twoLevel/ARG.cluster-x1";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer cluster1K file " << fileName << endl;
		exit(1);
	}
	index = 0;
	while(getline(file, s)){
		splitString(s, " ", &spt);
		int centerRotID = atoi(spt[0].c_str());
		for(int i=0;i<spt.size();i++){
			int rotID = atoi(spt[i].c_str());
			ARGIndex[rotID] = index;
			ARGCenter[rotID] = centerRotID;
		}
		index ++;
	}
	file.close();

	/*
	 * GLU, GLN, ARG, LYS, MET local energy change
	 */
	fileName = path + "scRotamer/rot/den/GLU.den";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer local density file " << fileName << endl;
		exit(1);
	}
	double GLUDensity[10000];
	index = 0;
	while(getline(file, s)){
		GLUDensity[index] = log(atof(s.c_str()));
		index++;
	}
	for(int i=0;i<10000;i++){
		GLUELocal[i] = GLUDensity[GLUCenter[i]] - GLUDensity[i];
	}
	file.close();

	fileName = path + "scRotamer/rot/den/GLN.den";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer local density file " << fileName << endl;
		exit(1);
	}
	double GLNDensity[10000];
	index = 0;
	while(getline(file, s)){
		GLNDensity[index] = log(atof(s.c_str()));
		index++;
	}
	for(int i=0;i<10000;i++){
		GLNELocal[i] = GLNDensity[GLNCenter[i]] - GLNDensity[i];
	}
	file.close();

	fileName = path + "scRotamer/rot/den/LYS.den";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer local density file " << fileName << endl;
		exit(1);
	}
	double LYSDensity[10000];
	index = 0;
	while(getline(file, s)){
		LYSDensity[index] = log(atof(s.c_str()));
		index++;
	}
	for(int i=0;i<10000;i++){
		LYSELocal[i] = LYSDensity[LYSCenter[i]] - LYSDensity[i];
	}
	file.close();

	fileName = path + "scRotamer/rot/den/MET.den";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer local density file " << fileName << endl;
		exit(1);
	}
	double METDensity[10000];
	index = 0;
	while(getline(file, s)){
		METDensity[index] = log(atof(s.c_str()));
		index++;
	}
	for(int i=0;i<10000;i++){
		METELocal[i] = METDensity[METCenter[i]] - METDensity[i];
	}
	file.close();

	fileName = path + "scRotamer/rot/den/ARG.den";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open sidechain rotamer local density file " << fileName << endl;
		exit(1);
	}
	double ARGDensity[20000];
	index = 0;
	while(getline(file, s)){
		ARGDensity[index] = log(atof(s.c_str()));
		index++;
	}
	for(int i=0;i<20000;i++){
		ARGELocal[i] = ARGDensity[ARGCenter[i]] - ARGDensity[i];
	}
	file.close();


	/*
	 * GLU, GLN, ARG, LYS, MET cluster: coverage*4
	 */

	for(int i=0;i<20;i++){
		for(int j=0;j<1000;j++){
			for(int k=0;k<210;k++){
				rotCluster[i][j][k] = 0;
			}
		}
	}

	vector<string> triList;
	triList.push_back("GLU");
	triList.push_back("GLN");
	triList.push_back("LYS");
	triList.push_back("MET");
	triList.push_back("ARG");
	for(int i=0;i<triList.size();i++){
		string tri = triList[i];
		fileName = path + "scRotamer/rot/twoLevel/"+tri+".cluster-x4";
		file.open(fileName.c_str(), ios::in);
		if(! file.is_open()) {
			cout << "fail to open sidechain rotamer cluster file " << fileName << endl;
			exit(1);
		}
		index = 0;
		while(getline(file, s)){
			splitString(s, " ", &spt);
			int centerRotID = atoi(spt[0].c_str());
			rotCluster[rn.triToInt(tri)][index][0] = spt.size();
			for(int i=0;i<spt.size();i++){
				int rotID = atoi(spt[i].c_str());
				rotCluster[rn.triToInt(tri)][index][i+1] = rotID;
			}
			index ++;
		}
		file.close();

        fileName = path + "scRotamer/rot/twoLevel/"+tri+".cluster-x1";
        file.open(fileName.c_str(), ios::in);
        if(! file.is_open()) {
                cout << "fail to open sidechain rotamer cluster file " << fileName << endl;
                exit(1);
        }
        index = 0;
        while(getline(file, s)){
                splitString(s, " ", &spt);
                int centerRotID = atoi(spt[0].c_str());
                rotClusterUnique[rn.triToInt(tri)][index][0] = spt.size();
                for(int i=0;i<spt.size();i++){
                        int rotID = atoi(spt[i].c_str());
                        rotClusterUnique[rn.triToInt(tri)][index][i+1] = rotID;
                }
                index ++;
        }
        file.close();
	}
}

int ResScRotamerLib::getRotamerID(ResScRotamer* rot){
	double minD = 99.9;
	double d;
	int minIndex = 0;
	for(int i=0;i<rotList[rot->aaType].size();i++){
		d = rotList[rot->aaType][i]->distanceTo(rot);
		if(d < minD){
			minD = d;
			minIndex = i;
		}
	}
	return minIndex;
}

int ResScRotamerLib::getAminoAcidRotamerClusterNum(int aa){
	//if(aa == 2 || aa == 11 || aa == 7 || aa == 9) return 50;
	//else if(aa == 6 || aa == 4 || aa == 18 || aa == 19) return 100;
	//else if(aa == 3 || aa == 13 || aa == 8 || aa == 10 || aa == 14) return 200;
	if(aa == 3 || aa == 13 || aa == 8 || aa == 10 || aa == 14) return 1000;
	else return 0;
}

void ResScRotamerLib::checkResScRotamerPolarGroups() {
	for(int i=0;i<20;i++){
		ResScRotamer* rot = rotList[i][0];
		int polarAtomNum = rot->polarAtomNum;
		for(int j=0;j<polarAtomNum;j++) {
			int uniqueID = rot->uniqueIDs[rot->polarAtomIndex[j]];

		}
	}
}

ResScRotamerLib::~ResScRotamerLib() {
	// TODO Auto-generated destructor stub
	delete this->atLib;
	for(int i=0;i<20;i++){
		for(int j=0;j<rotList[i].size();j++){
			delete rotList[i][j];
		}
	}
}

}
