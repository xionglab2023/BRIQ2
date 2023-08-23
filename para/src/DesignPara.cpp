/*
 * DesignPara.cpp
 *
 *  Created on: 2022Äê3ÔÂ10ÈÕ
 *      Author: pengx
 */

#include "para/DesignPara.h"

namespace NSPpara {

DesignPara::DesignPara(const string& file) {
	// TODO Auto-generated constructor stub
	this->atLib = new AtomLib();
	ifstream f;
	f.open(file, ios::in);
	if(!f.is_open()){
		cout << "fail to open file: " << endl;
		exit(1);
	}

	this->shiftL = 0.1;
	this->shiftS = 0.1;

	string s;
	int id = 0;
	for(int i=0;i<60;i++){
		ref[i] = 0.0;
	}

	while(getline(f, s)){
		if(id < 5){
			vdwSepWeight[id] = atof(s.c_str());
		}
		else if(id == 5){
			curve = s;
		}
		else if(id == 6){
			vdwRange = atof(s.c_str());
		}
		else if(id == 7)
			wtRot = atof(s.c_str());
		else if(id == 8) {
			lamdaS = atof(s.c_str());
		}
		else if(id == 9) {
			lamdaL = atof(s.c_str());
		}
		else if(id == 10) {
			this->wdHb = atof(s.c_str());
		}
		else if(id == 11) {
			this->polarKS = atof(s.c_str());
		}
		else if(id == 12) {
			this->polarKL = atof(s.c_str());
		}
		else if(id == 13) {
			this->polarCutS = atof(s.c_str());
		}
		else if(id == 14){
			this->polarCutL = atof(s.c_str());
		}
		else if(id == 15){
			wd0 = atof(s.c_str());
		}
		else if(id < 76) {
			vdwWdRescale[id-16] = atof(s.c_str());
		}
		else if(id < 80) {
			polarWdRescale[id-76] = atof(s.c_str());
		}
		else if(id < 150) {
			polarWd[id-80] = -1.0* atof(s.c_str());
		}
		else if(id < 210) {
			ref[id-150] = atof(s.c_str());
		}
		id++;
	}
	f.close();


	this->wS1 = 1.0;

	this->wS245[ 0] = 0.7448;
	this->wS245[ 1] = 0.5725;
	this->wS245[ 2] = 0.7695;
	this->wS245[ 3] = 0.8173;
	this->wS245[ 4] = 0.7281;
	this->wS245[ 5] = 0.1774;
	this->wS245[ 6] = 0.5349;
	this->wS245[ 7] = 0.2435;
	this->wS245[ 8] = 0.1461;
	this->wS245[ 9] = 0.6280;
	this->wS245[10] = 0.5951;
	this->wS245[11] = 0.4741;
	this->wS245[12] = 0.7619;
	this->wS245[13] = 0.6206;
	this->wS245[14] = 0.8467;
	this->wS245[15] = 0.7039;
	this->wS245[16] = 0.4684;
	this->wS245[17] = 0.1720;
	this->wS245[18] = 0.3304;
	this->wS245[19] = 0.5774;
	this->wS245[20] = 0.5064;
	this->wS245[21] = 0.6610;
	this->wS245[22] = 0.2828;
	this->wS245[23] = 0.5887;
	this->wS245[24] = 0.5354;
	this->wS245[25] = 0.4182;
	this->wS245[26] = 0.6699;
	this->wS245[27] = 0.5874;
	this->wS245[28] = 0.4612;
	this->wS245[29] = 0.7305;
	this->wS245[30] = 0.3695;
	this->wS245[31] = 0.5406;
	this->wS245[32] = 0.5982;
	this->wS245[33] = 0.6971;
	this->wS245[34] = 0.9017;
	this->wS245[35] = 0.5244;
	this->wS245[36] = 0.4347;
	this->wS245[37] = 0.6456;
	this->wS245[38] = 0.4433;
	this->wS245[39] = 0.7392;
	this->wS245[40] = 0.6689;
	this->wS245[41] = 0.7525;
	this->wS245[42] = 0.5548;
	this->wS245[43] = 0.5492;
	this->wS245[44] = 0.9816;

	this->dsPara[0] = 90.0; //disulfide bond length
	this->dsPara[1]	= 0.02; //disulfide bond angle
	this->dsPara[2] = 1.0; //disulfide bond dihedral angle
	this->dsPara[3] = -20.0; //well depth

	string path = NSPdataio::datapath();
	string nbDistanceFile = path+"/atomLib/nb.dist";

	for(int i=0;i<167;i++){
		for(int j=0;j<167;j++){
			this->vdwNbD0[i][j] = atLib->getAtomProperty(i)->vdwRadius + atLib->getAtomProperty(j)->vdwRadius;
		}
	}

	ResName rn;
	ifstream nbfile;
	nbfile.open(nbDistanceFile.c_str(), ios::in);
	vector<string> spt;
	while(getline(nbfile, s)){
		NSPtools::splitString(s, " ", &spt);
		int idA = atLib->uniqueNameToUniqueID[spt[0]];
		for(int i=0;i<20;i++){
			int idB = atLib->uniqueNameToID(rn.intToTri(i)+"-"+spt[1]);
			this->vdwNbD0[idA][idB] = atof(spt[2+i].c_str());
			this->vdwNbD0[idB][idA] = atof(spt[2+i].c_str());
		}
	}
	nbfile.close();

	updateVdwWellDepthDesign();
}

void DesignPara::updateVdwWellDepthDesign() {
	map<int, int> pairIndexToVdwParaIndex;
	int index = 0;
	for(int i=0;i<5;i++){
		for(int j=i;j<5;j++){
			pairIndexToVdwParaIndex[i*5+j] = index;
			pairIndexToVdwParaIndex[j*5+i] = index;
			index++;
		}
	}

	/*
	for(int i=0;i<167;i++){
		int typeA = atLib->uniqueIDToVdwAtomID[i];
		int nbA = atLib->apList[i]->connectNum;
		for(int j=0;j<167;j++){
			int typeB = atLib->uniqueIDToVdwAtomID[j];
			int paraIndex = pairIndexToVdwParaIndex[typeA*5+typeB];
			int nbB = atLib->apList[j]->connectNum;

			this->vdwRescale1[i][j] = this->vdwWdRescale[paraIndex];
			this->vdwRescale2[i][j] = this->vdwWdRescale[paraIndex + 15];
			this->vdwRescale3[i][j] = this->vdwWdRescale[paraIndex + 30];
			this->vdwRescale4[i][j] = this->vdwWdRescale[paraIndex + 45];

		}
	}
	*/
}

DesignPara::~DesignPara() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPpara */
