/*
 * ProParameter.cpp
 *
 */


#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "tools/StringTool.h"
#include "para/ProParameter.h"

namespace NSPpara {

ProParameter::ProParameter() {
	// TODO Auto-generated constructor stub

	this->curve = "6";

	this->vdwRange = 6.55;

	this->wSurf = 0.05;
	this->wCore = 0.72;
	this->saiD0 = 0.65;
	this->saiSigma = 0.07;

	this->wS1 = 1.0;
	this->wS2 = 1.0;
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


	wtRot = 1.0;

	for(int i=0;i<20;i++) {
		this->ref[i] = 0.0;
	}

	double vol[14];
	int polarType[14];
	double polarVdwRescale[6];
	polarVdwRescale[0] = 1.0;
	polarVdwRescale[1] = 1.0;
	polarVdwRescale[2] = 0.5;
	polarVdwRescale[3] = 1.0;
	polarVdwRescale[4] = 0.5;
	polarVdwRescale[5] = 0.0;

	polarType[0] = 0; //CA
	polarType[1] = 0; //CT1
	polarType[2] = 0; //CT2
	polarType[3] = 0; //CT3
	polarType[4] = 1; //CR0
	polarType[5] = 1; //CR1
	polarType[6] = 1; //CP0
	polarType[7] = 2; //O1
	polarType[8] = 2; //O2
	polarType[9] = 2; //O3
	polarType[10] = 2; //N1
	polarType[11] = 2; //N2
	polarType[12] = 2; //N3
	polarType[13] = 0; //S

	double pairPolarVdwRescale[105];
	int volIndex = 0;
	double wt;
	for(int i=0;i<14;i++){
		int typeA = polarType[i];
		for(int j=i;j<14;j++){
			int typeB = polarType[j];
			if(typeA == 0 && typeB == 0)
				wt = polarVdwRescale[0];
			else if(typeA == 0 && typeB == 1)
				wt = polarVdwRescale[1];
			else if(typeA == 0 && typeB == 2)
				wt = polarVdwRescale[2];
			else if(typeA == 1 && typeB == 0)
				wt = polarVdwRescale[1];
			else if(typeA == 1 && typeB == 1)
				wt = polarVdwRescale[3];
			else if(typeA == 1 && typeB == 2)
				wt = polarVdwRescale[4];
			else if(typeA == 2 && typeB == 0)
				wt = polarVdwRescale[2];
			else if(typeA == 2 && typeB == 1)
				wt = polarVdwRescale[4];
			else if(typeA == 2 && typeB == 2)
				wt = polarVdwRescale[5];

			pairPolarVdwRescale[volIndex] = wt;
			volIndex++;
		}
	}

	this->vdwNbRescale[0] = 0.6667;
	this->vdwNbRescale[1] = 0.333;
	this->vdwNbRescale[2] = 0.25;
	this->vdwNbRescale[3] = 0.2;
	this->vdwNbRescale[4] = 0.1667;


	for(int i=0;i<105;i++){
		paraVdw[i*5] = 0.0;
		paraVdw[i*5+1] = 2.5;
		paraVdw[i*5+2] = 0.0;
		paraVdw[i*5+3] = 2.5;
		paraVdw[i*5+4] = -1.3 * pairPolarVdwRescale[i];
	}

	vdwSepWeight[0] = 0.5;
	vdwSepWeight[1] = 1.0;
	vdwSepWeight[2] = 1.0;
	vdwSepWeight[3] = 1.0;
	vdwSepWeight[4] = 1.0;

	polarSepWeight[0] = 0.0;
	polarSepWeight[1] = 1.0;
	polarSepWeight[2] = 1.0;
	polarSepWeight[3] = 1.0;
	polarSepWeight[4] = 1.0;


	char xx[20];
	this->atLib = new AtomLib();

	updateVdwWellDepthDesign();

	for(int i=0;i<238;i++){
		this->paraPolar[i][0] = 0.0; //polar well depth
		this->paraPolar[i][1] = 5.0; //density range
		this->paraPolar[i][2] = 0.0; //k2

		this->paraPolar[i][3] = 0.0; //polar well depth
		this->paraPolar[i][4] = 5.0; //k
		this->paraPolar[i][5] = 0.0; //k2
	}

	for(int i=0;i<14;i++){
		this->polarDamping[i] = 3.7;
	}

	for(int i=14;i<28;i++){
		this->polarDamping[i] = 0.07;
	}

	this->dsPara[0] = 90.0; //disulfide bond length
	this->dsPara[1]	= 0.02; //disulfide bond angle
	this->dsPara[2] = 1.0; //disulfide bond dihedral angle
	this->dsPara[3] = -31.0; //well depth

	string path = NSPdataio::datapath();
	string nbDistanceFile = path+"/atomLib/nb.dist";

	for(int i=0;i<167;i++){
		for(int j=0;j<167;j++){
			this->vdwNbD0[i][j] = atLib->getAtomProperty(i)->vdwRadius + atLib->getAtomProperty(j)->vdwRadius;
		}
	}

	ResName rn;
	ifstream file;
	file.open(nbDistanceFile.c_str(), ios::in);
	string s;
	vector<string> spt;
	while(getline(file, s)){
		splitString(s, " ", &spt);
		int idA = atLib->uniqueNameToUniqueID[spt[0]];
		for(int i=0;i<20;i++){
			int idB = atLib->uniqueNameToID(rn.intToTri(i)+"-"+spt[1]);
			this->vdwNbD0[idA][idB] = atof(spt[2+i].c_str());
			this->vdwNbD0[idB][idA] = atof(spt[2+i].c_str());
		}
	}
	file.close();


}

/*
ProParameter::ProParameter(int type, const string& file, const string& filePolar){

	 //* type = 0: parameter for packing
	 //* type = 1: parameter for design
	 //* vdw well depth is volumn dependent

	this->atLib = new ProteinAtomLib();
	this->dsPara[0] = 90.0; //disulfide bond length
	this->dsPara[1]	= 0.02; //disulfide bond angle
	this->dsPara[2] = 1.0; //disulfide bond dihedral angle
	this->dsPara[3] = -31.0; //well depth

	double polarVdwRescale[6];
	for(int i=0;i<6;i++){
		polarVdwRescale[i] = 1.0;
	}

	this->wS245[ 0] = 0.2338;
	this->wS245[ 1] = 0.3320;
	this->wS245[ 2] = 0.2793;
	this->wS245[ 3] = 0.4180;
	this->wS245[ 4] = 0.7895;
	this->wS245[ 5] = 0.3504;
	this->wS245[ 6] = 0.3554;
	this->wS245[ 7] = 0.5124;
	this->wS245[ 8] = 0.5584;
	this->wS245[ 9] = 0.8110;
	this->wS245[10] = 0.2890;
	this->wS245[11] = 0.3112;
	this->wS245[12] = 0.2768;
	this->wS245[13] = 0.4983;
	this->wS245[14] = 0.8199;
	this->wS245[15] = 0.3104;
	this->wS245[16] = 0.6615;
	this->wS245[17] = 0.5480;
	this->wS245[18] = 0.4350;
	this->wS245[19] = 0.6700;
	this->wS245[20] = 0.1558;
	this->wS245[21] = 0.2162;
	this->wS245[22] = 0.4735;
	this->wS245[23] = 1.0074;
	this->wS245[24] = 0.5618;
	this->wS245[25] = 0.2341;
	this->wS245[26] = 0.4209;
	this->wS245[27] = 0.6252;
	this->wS245[28] = 0.4096;
	this->wS245[29] = 0.6972;
	this->wS245[30] = 0.2002;
	this->wS245[31] = 0.4907;
	this->wS245[32] = 0.6281;
	this->wS245[33] = 0.5749;
	this->wS245[34] = 0.8163;
	this->wS245[35] = 0.4122;
	this->wS245[36] = 0.3237;
	this->wS245[37] = 0.5074;
	this->wS245[38] = 0.4578;
	this->wS245[39] = 0.7095;
	this->wS245[40] = 0.2285;
	this->wS245[41] = 0.4443;
	this->wS245[42] = 0.3358;
	this->wS245[43] = 0.3803;
	this->wS245[44] = 0.9797;

	if(type == 1) {
		ifstream f;
		f.open(file, ios::in);
		if(!f.is_open()){
			cout << "fail to open file: " << file << endl;
			exit(1);
		}

		string s;
		int id = 0;
		double kS = 0;
		double shiftS = 0;
		double kL = 0;
		double shiftL = 0;
		double kVol = 2.0;

		while(getline(f, s)){
			if(id < 5){
				vdwSepWeight[id] = atof(s.c_str());
			}
			else if(id < 10){
				polarSepWeight[id-5] = atof(s.c_str());
			}
			else if(id == 10){
				this->curve = s;
			}
			else if(id == 11)
				kS = atof(s.c_str());
			else if(id == 12)
				shiftS = atof(s.c_str());
			else if(id == 13)
				kL = atof(s.c_str());
			else if(id == 14)
				shiftL = atof(s.c_str());
			else if(id == 15){
				vdwRange = atof(s.c_str());
			}
			else if(id == 16){
				saiD0 = atof(s.c_str());
			}
			else if(id == 17){
				saiSigma = atof(s.c_str());
			}
			else if(id == 18) {
				wSurf = atof(s.c_str());
			}
			else if(id == 19) {
				double wr = atof(s.c_str());
				for(int i=0;i<20;i++){
					wtRot[i] = wr;
				}
			}
			else if(id == 20) {
				wS1 = atof(s.c_str());
			}
			else if(id == 21) {
				wS2 = atof(s.c_str());
			}
			else if(id == 22) {
				kVol = atof(s.c_str());
			}
			else if(id >= 23 && id < 29) {
				polarVdwRescale[id-23] = atof(s.c_str());
			}
			else if(id >= 29 && id < 34) {
				wS2Sep[id-29] = atof(s.c_str());
			}
			id++;
		}

		if(id != 34){
			cout << "parameter file error: " << id << endl;
			exit(1);
		}
		f.close();

		double vol[14];
		int polarType[14];
		vol[0] = 0.38;   polarType[0] = 0; //CA
		vol[1] = 0.38;   polarType[1] = 0; //CT1
		vol[2] = 0.63;   polarType[2] = 0; //CT2
		vol[3] = 1.00;   polarType[3] = 0; //CT3
		vol[4] = 0.22;   polarType[4] = 1; //CR0
		vol[5] = 0.55;   polarType[5] = 1; //CR1
		vol[6] = 0.22;   polarType[6] = 1; //CP0
		vol[7] = 0.45;   polarType[7] = 2; //O1
		vol[8] = 0.52;   polarType[8] = 2; //O2
		vol[9] = 0.52;   polarType[9] = 2; //O3
		vol[10] = 0.33;  polarType[10] = 2; //N1
		vol[11] = 1.05;  polarType[11] = 2; //N2
		vol[12] = 0.70;  polarType[12] = 2; //N3
		vol[13] = 0.80;  polarType[13] = 0; //S

		double pairVol[105];
		int volIndex = 0;
		double wt;
		for(int i=0;i<14;i++){
			int typeA = polarType[i];
			for(int j=i;j<14;j++){
				int typeB = polarType[j];
				if(typeA == 0 && typeB == 0)
					wt = polarVdwRescale[0];
				else if(typeA == 0 && typeB == 1)
					wt = polarVdwRescale[1];
				else if(typeA == 0 && typeB == 2)
					wt = polarVdwRescale[2];
				else if(typeA == 1 && typeB == 0)
					wt = polarVdwRescale[1];
				else if(typeA == 1 && typeB == 1)
					wt = polarVdwRescale[3];
				else if(typeA == 1 && typeB == 2)
					wt = polarVdwRescale[4];
				else if(typeA == 2 && typeB == 0)
					wt = polarVdwRescale[2];
				else if(typeA == 2 && typeB == 1)
					wt = polarVdwRescale[4];
				else if(typeA == 2 && typeB == 2)
					wt = polarVdwRescale[5];

				pairVol[volIndex] = wt*pow(vol[i]*vol[j], kVol);
				volIndex++;
			}
		}
		if(volIndex != 105){
			cout << "vol index error" << volIndex << endl;
			exit(1);
		}

		for(int i=0;i<105;i++){
			paraVdw[i*5] = shiftS;
			paraVdw[i*5+1] = kS;
			paraVdw[i*5+2] = shiftL;
			paraVdw[i*5+3] = kL;
			paraVdw[i*5+4] = -2.0*pairVol[i];
		}

		updateVdwWellDepth();

		f.open(filePolar, ios::in);
		if(!f.is_open()){
			cout << "fail to open file: " << endl;
			exit(1);
		}

		for(int i=0;i<14;i++){
			this->polarDamping[i] = 3.7;
		}

		polarDamping[2] = 3.6;
		polarDamping[3] = 3.6;
		polarDamping[4] = 3.6;
		polarDamping[5] = 3.6;
		polarDamping[6] = 3.9;
		polarDamping[7] = 3.6;
		polarDamping[10] = 3.6;
		polarDamping[11] = 3.9;
		polarDamping[12] = 3.9;
		polarDamping[13] = 3.9;

		for(int i=14;i<28;i++){
			this->polarDamping[i] = 0.07;
		}


		vector<string> spt;
		id = 0;
		while(getline(f, s)){
			if(id < 1428){
				this->paraPolar[id/6][id%6] = atof(s.c_str());
			}
			id++;
		}

	}
	else {
		cout << "invalid para type" << endl;
		exit(1);
	}
}
*/

ProParameter::ProParameter(bool designTag, const string& file, const string& filePolar){
	/*
	 * vdw well depth is neighbor number dependent
	 */
	this->atLib = new AtomLib();
	this->dsPara[0] = 90.0; //disulfide bond length
	this->dsPara[1]	= 0.02; //disulfide bond angle
	this->dsPara[2] = 1.0; //disulfide bond dihedral angle
	this->dsPara[3] = -31.0; //well depth

	double polarVdwRescale[6];
	for(int i=0;i<6;i++){
		polarVdwRescale[i] = 1.0;
	}

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

	if(designTag == true) {
		ifstream f;
		f.open(file, ios::in);
		if(!f.is_open()){
			cout << "fail to open file: " << file << endl;
			exit(1);
		}

		string s;
		int id = 0;
		double kS = 0;
		double shiftS = 0;
		double kL = 0;
		double shiftL = 0;

		while(getline(f, s)){
			if(id < 5){
				vdwSepWeight[id] = atof(s.c_str());
			}
			else if(id < 10){
				polarSepWeight[id-5] = atof(s.c_str());
			}
			else if(id == 10){
				this->curve = s;
			}
			else if(id == 11)
				kS = atof(s.c_str());
			else if(id == 12)
				shiftS = atof(s.c_str());
			else if(id == 13)
				kL = atof(s.c_str());
			else if(id == 14)
				shiftL = atof(s.c_str());
			else if(id == 15){
				vdwRange = atof(s.c_str());
			}
			else if(id == 16){
				saiD0 = atof(s.c_str());
			}
			else if(id == 17){
				saiSigma = atof(s.c_str());
			}
			else if(id == 18) {
				wSurf = atof(s.c_str());
			}
			else if(id == 19) {
				double wr = atof(s.c_str());
				wtRot = wr;
			}
			else if(id == 20) {
				wS1 = atof(s.c_str());
			}
			else if(id == 21) {
				wS2 = atof(s.c_str());
			}
			else if(id >= 22 && id < 27) {
				vdwNbRescale[id-22] = atof(s.c_str());
			}
			else if(id >= 27 && id < 33) {
				polarVdwRescale[id-27] = atof(s.c_str());
			}
			else if(id >=33 && id < 53) {
				ref[id-33] = atof(s.c_str());
			}
			id++;
		}

		if(id != 53){
			cout << "parameter file error: " << id << endl;
			exit(1);
		}
		f.close();

		double vol[14];
		int polarType[14];
		vol[0] = 0.38;   polarType[0] = 0; //CA
		vol[1] = 0.38;   polarType[1] = 0; //CT1
		vol[2] = 0.63;   polarType[2] = 0; //CT2
		vol[3] = 1.00;   polarType[3] = 0; //CT3
		vol[4] = 0.22;   polarType[4] = 1; //CR0
		vol[5] = 0.55;   polarType[5] = 1; //CR1
		vol[6] = 0.22;   polarType[6] = 1; //CP0
		vol[7] = 0.45;   polarType[7] = 2; //O1
		vol[8] = 0.52;   polarType[8] = 2; //O2
		vol[9] = 0.52;   polarType[9] = 2; //O3
		vol[10] = 0.33;  polarType[10] = 2; //N1
		vol[11] = 1.05;  polarType[11] = 2; //N2
		vol[12] = 0.70;  polarType[12] = 2; //N3
		vol[13] = 0.80;  polarType[13] = 0; //S

		double pairPolarVdwRescale[105];
		int volIndex = 0;
		double wt;
		for(int i=0;i<14;i++){
			int typeA = polarType[i];
			for(int j=i;j<14;j++){
				int typeB = polarType[j];
				if(typeA == 0 && typeB == 0)
					wt = polarVdwRescale[0];
				else if(typeA == 0 && typeB == 1)
					wt = polarVdwRescale[1];
				else if(typeA == 0 && typeB == 2)
					wt = polarVdwRescale[2];
				else if(typeA == 1 && typeB == 0)
					wt = polarVdwRescale[1];
				else if(typeA == 1 && typeB == 1)
					wt = polarVdwRescale[3];
				else if(typeA == 1 && typeB == 2)
					wt = polarVdwRescale[4];
				else if(typeA == 2 && typeB == 0)
					wt = polarVdwRescale[2];
				else if(typeA == 2 && typeB == 1)
					wt = polarVdwRescale[4];
				else if(typeA == 2 && typeB == 2)
					wt = polarVdwRescale[5];

				pairPolarVdwRescale[volIndex] = wt;
				volIndex++;
			}
		}
		if(volIndex != 105){
			cout << "vol index error" << volIndex << endl;
			exit(1);
		}

		for(int i=0;i<105;i++){
			paraVdw[i*5] = shiftS;
			paraVdw[i*5+1] = kS;
			paraVdw[i*5+2] = shiftL;
			paraVdw[i*5+3] = kL;
			paraVdw[i*5+4] = -2.0*pairPolarVdwRescale[i];
		}

		updateVdwWellDepthDesign();

		f.open(filePolar, ios::in);
		if(!f.is_open()){
			cout << "fail to open file: " << endl;
			exit(1);
		}

		for(int i=0;i<14;i++){
			this->polarDamping[i] = 3.7;
		}

		polarDamping[2] = 3.6;
		polarDamping[3] = 3.6;
		polarDamping[4] = 3.6;
		polarDamping[5] = 3.6;
		polarDamping[6] = 3.9;
		polarDamping[7] = 3.6;
		polarDamping[10] = 3.6;
		polarDamping[11] = 3.9;
		polarDamping[12] = 3.9;
		polarDamping[13] = 3.9;

		for(int i=14;i<28;i++){
			this->polarDamping[i] = 0.07;
		}

		vector<string> spt;
		id = 0;
		while(getline(f, s)){
			if(id < 1428){
				this->paraPolar[id/6][id%6] = atof(s.c_str());
			}
			id++;
		}
	}
	else {
		cout << "invalid para type" << endl;
		exit(1);
	}

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
	string s;
	vector<string> spt;
	while(getline(nbfile, s)){
		splitString(s, " ", &spt);
		int idA = atLib->uniqueNameToUniqueID[spt[0]];
		for(int i=0;i<20;i++){
			int idB = atLib->uniqueNameToID(rn.intToTri(i)+"-"+spt[1]);
			this->vdwNbD0[idA][idB] = atof(spt[2+i].c_str());
			this->vdwNbD0[idB][idA] = atof(spt[2+i].c_str());
		}
	}
	nbfile.close();
}



ProParameter::ProParameter(const string& file){

	/*
	 * vdw para: 0~5
	 * wtRot: 6
	 * curve type: 7
	 * volFactor: 8
	 */

	this->curve = "3";

	this->dsPara[0] = 90.0; //disulfide bond length
	this->dsPara[1]	= 0.02; //disulfide bond angle
	this->dsPara[2] = 1.0; //disulfide bond dihedral angle
	this->dsPara[3] = -31.0; //well depth
	wtRot = 2.2;

	polarSepWeight[0] = 0.0;
	polarSepWeight[1] = 1.0;
	polarSepWeight[2] = 1.0;
	polarSepWeight[3] = 1.0;
	polarSepWeight[4] = 1.0;

	ifstream f;
	f.open(file, ios::in);
	if(!f.is_open()){
		cout << "fail to open file: " << endl;
		exit(1);
	}

	string s;
	int id = 0;
	while(getline(f, s)){
		if(id < 5){
			vdwSepWeight[id] = atof(s.c_str());
		}
		else if(id < 530) {
			paraVdw[id-5] = atof(s.c_str());
		}
		id++;
	}
	f.close();

	char xx[20];
	this->atLib = new AtomLib();

	updateVdwWellDepth();

	for(int i=0;i<238;i++){
		this->paraPolar[i][0] = -1.0; //polar well depth
		this->paraPolar[i][1] = 5.0; //density range
		this->paraPolar[i][2] = 0.2; //k2

		this->paraPolar[i][3] = -1.0; //polar well depth
		this->paraPolar[i][4] = 5.0; //k
		this->paraPolar[i][5] = 0.2; //k2
	}

	for(int i=0;i<14;i++){
		this->polarDamping[i] = 3.7;
	}

	for(int i=14;i<28;i++){
		this->polarDamping[i] = 0.07;
	}

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
		splitString(s, " ", &spt);
		int idA = atLib->uniqueNameToUniqueID[spt[0]];
		for(int i=0;i<20;i++){
			int idB = atLib->uniqueNameToID(rn.intToTri(i)+"-"+spt[1]);
			this->vdwNbD0[idA][idB] = atof(spt[2+i].c_str());
			this->vdwNbD0[idB][idA] = atof(spt[2+i].c_str());
		}
	}
	nbfile.close();
}

ProParameter::ProParameter(const string& fileVdw, const string& filePolar){

	/*
	 * vdw para: 0~5
	 * wtRot: 6
	 * curve type: 7
	 * volFactor: 8
	 */
	this->curve = "3";

	ifstream f;
	f.open(fileVdw, ios::in);
	if(!f.is_open()){
		cout << "fail to open file: " << fileVdw << endl;
		exit(1);
	}

	wtRot = 2.2;

	polarSepWeight[0] = 0.0;
	polarSepWeight[1] = 1.0;
	polarSepWeight[2] = 1.0;
	polarSepWeight[3] = 1.0;
	polarSepWeight[4] = 1.0;

	string s;
	int id = 0;
	while(getline(f, s)){
		if(id < 5){
			vdwSepWeight[id] = atof(s.c_str());
		}
		else if(id < 530) {
			paraVdw[id-5] = atof(s.c_str());
		}
		id++;
	}
	f.close();


	char xx[20];


	this->atLib = new AtomLib();
	updateVdwWellDepth();


	f.open(filePolar, ios::in);
	if(!f.is_open()){
		cout << "fail to open file: " << endl;
		exit(1);
	}

	for(int i=0;i<14;i++){
		this->polarDamping[i] = 3.7;
	}

	polarDamping[2] = 3.6;
	polarDamping[3] = 3.6;
	polarDamping[4] = 3.6;
	polarDamping[5] = 3.6;
	polarDamping[6] = 3.9;
	polarDamping[7] = 3.6;
	polarDamping[10] = 3.6;
	polarDamping[11] = 3.9;
	polarDamping[12] = 3.9;
	polarDamping[13] = 3.9;

	for(int i=14;i<28;i++){
		this->polarDamping[i] = 0.07;
	}


	vector<string> spt;
	id = 0;
	while(getline(f, s)){
		if(id < 1428){
			this->paraPolar[id/6][id%6] = atof(s.c_str());
		}
		id++;
	}

	this->dsPara[0] = 90.0; //disulfide bond length
	this->dsPara[1]	= 0.02; //disulfide bond angle
	this->dsPara[2] = 1.0; //disulfide bond dihedral angle
	this->dsPara[3] = -31.0; //well depth


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
	while(getline(nbfile, s)){
		splitString(s, " ", &spt);
		int idA = atLib->uniqueNameToUniqueID[spt[0]];
		for(int i=0;i<20;i++){
			int idB = atLib->uniqueNameToID(rn.intToTri(i)+"-"+spt[1]);
			this->vdwNbD0[idA][idB] = atof(spt[2+i].c_str());
			this->vdwNbD0[idB][idA] = atof(spt[2+i].c_str());
		}
	}
	nbfile.close();
}


void ProParameter::updateVdwWellDepthDesign() {
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

			this->vdwShiftS[i][j] = this->paraVdw[paraIndex*5];
			this->vdwLamdaS[i][j] = this->paraVdw[paraIndex*5+1];
			this->vdwShiftL[i][j] = this->paraVdw[paraIndex*5+2];
			this->vdwLamdaL[i][j] = this->paraVdw[paraIndex*5+3];
			this->vdwWd[i][j] = this->paraVdw[paraIndex*5+4] * this->vdwNbRescale[nbA+nbB-2];
		}
	}
	*/
}


void ProParameter::updateVdwWellDepth(){

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
		for(int j=0;j<167;j++){
			int typeB = atLib->uniqueIDToVdwAtomID[j];
			int paraIndex = pairIndexToVdwParaIndex[typeA*14+typeB];

			this->vdwShiftS[i][j] = this->paraVdw[paraIndex*5];
			this->vdwLamdaS[i][j] = this->paraVdw[paraIndex*5+1];
			this->vdwShiftL[i][j] = this->paraVdw[paraIndex*5+2];
			this->vdwLamdaL[i][j] = this->paraVdw[paraIndex*5+3];
			this->vdwWd[i][j] = this->paraVdw[paraIndex*5+4];
		}
	}
	*/
}

ProParameter::~ProParameter() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPmodel */
