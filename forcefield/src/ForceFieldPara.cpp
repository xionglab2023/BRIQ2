/*
 * ForceFieldPara.cpp
 *
 *  Created on: 2023Äê4ÔÂ25ÈÕ
 *      Author: pengx
 */

#include "forcefield/ForceFieldPara.h"

namespace NSPforcefield {

ForceFieldPara::ForceFieldPara() {
	// TODO Auto-generated constructor stub
	this->clashLamdaNb = 2.5;
	this->clashLamdaNnb = 2.5;

	this->hbLamda1 = 0.1;
	this->hbLamda2 = 0.15;
	this->wtHb = 1.0;

	this->rnaKAng = 0.1;
	this->rnaKBond = 7.0;

	this->wtPho = 0.7;
	this->wtRibose = 1.0;

	for(int i=0;i<6;i++){
		rnaDihedImpD1D2Shift[i] = 0.0;
		rnaDihedImpD4D5Shift[i] = 0.0;
		rnaDihedD2D3D4Shift[i] = 0.0;
	}

	for(int i=0;i<16;i++){
		rnaRiboseRotamerShift[i] = 0.0;
	}


	rnaRiboseRotamerShift[0] = 1.3;   rnaRiboseRotamerShift[1] = -1.8;
	rnaRiboseRotamerShift[4] = 1.4;   rnaRiboseRotamerShift[5] = -2.0;
	rnaRiboseRotamerShift[8] = 2.3;   rnaRiboseRotamerShift[9] = -1.5;
	rnaRiboseRotamerShift[12] = 2.4;  rnaRiboseRotamerShift[13] = -1.5;


	for(int i=0;i<16;i++){
		wtRiboseOxy[i] = 1.0;
	}


	this->wtBp1 = 1.0;
	this->wtBp2 = 1.0;
//	this->bbClash = 0.4;

	this->wtO4O2C2Nb = 0.7;
	this->wtO4O2C2Nnb = 0.7;

	this->wtClash = 1.0;

	this->T0 = 2.5;
	this->T1 = 0.1;
	this->anneal = 0.95;
	this->stepNum = 10;
	this->outFreq = 1000;

	this->lamdaClash = 2.8;

	this->initConnectWT = 0.4;
	this->initClashWT = 0.002;

	this->initShift = 0.6;
	this->dShift = 0.02;

	this->connectWTFactor = 1.05;
	this->clashWTFactor = 1.05;

	loopRiboConnectMove = true;
	ctRandMove = true;
	f3Move = true;
	singleBaseMove = true;
	reverseRotMove = true;

	this->spType = "sp2000";
	this->phoRep = 4.0;

}


ForceFieldPara::ForceFieldPara(const string& paraFile){
	this->clashLamdaNb = 2.5;
	this->clashLamdaNnb = 2.5;

	this->hbLamda1 = 0.1;
	this->hbLamda2 = 0.15;
	this->wtHb = 1.0;

	this->rnaKAng = 0.1;
	this->rnaKBond = 7.0;

	this->wtPho = 0.7;
	this->wtRibose = 1.0;

	for(int i=0;i<6;i++){
		rnaDihedImpD1D2Shift[i] = 0.0;
		rnaDihedImpD4D5Shift[i] = 0.0;
		rnaDihedD2D3D4Shift[i] = 0.0;
	}

	for(int i=0;i<16;i++){
		rnaRiboseRotamerShift[i] = 0.0;
	}

	rnaRiboseRotamerShift[0] = 1.3;   rnaRiboseRotamerShift[1] = -1.8;
	rnaRiboseRotamerShift[4] = 1.4;   rnaRiboseRotamerShift[5] = -2.0;
	rnaRiboseRotamerShift[8] = 2.3;   rnaRiboseRotamerShift[9] = -1.5;
	rnaRiboseRotamerShift[12] = 2.4;  rnaRiboseRotamerShift[13] = -1.5;


	for(int i=0;i<16;i++){
		wtRiboseOxy[i] = 1.0;
	}


	this->wtBp1 = 1.0;
	this->wtBp2 = 1.0;

	this->wtO4O2C2Nb = 0.7;
	this->wtO4O2C2Nnb = 0.7;

	this->wtClash = 1.0;

	this->T0 = 2.5;
	this->T1 = 0.1;
	this->anneal = 0.95;
	this->stepNum = 10;
	this->outFreq = 1000;

	this->lamdaClash = 2.8;

	this->initConnectWT = 0.4;
	this->initClashWT = 0.002;

	this->initShift = 0.6;
	this->dShift = 0.02;

	this->connectWTFactor = 1.05;
	this->clashWTFactor = 1.05;

	loopRiboConnectMove = true;
	ctRandMove = true;
	f3Move = true;
	singleBaseMove = true;
	reverseRotMove = true;

	this->spType = "sp2000";
	this->phoRep = 4.0;


	ifstream input;
	input.open(paraFile.c_str(), ios::in);
	if(!input.is_open()) {
		cout << "fail to open file: " << paraFile << endl;
		exit(1);
	}
	string s;
	vector<string> spt;
	vector<double> pList;
	while(getline(input,s)) {
		pList.push_back(atof(s.c_str()));
	}
	this->clashLamdaNb = pList[0];
	this->clashLamdaNnb = pList[1];
	this->hbLamda1 = pList[2];
	this->hbLamda2 = pList[3];
	this->wtHb = pList[4];
	this->rnaKBond = pList[5];
	this->rnaKAng = pList[6];
	this->wtPho = pList[7];
	//this->wtRibose = pList[8];
	for(int i=0;i<6;i++){
		rnaDihedImpD1D2Shift[i] = pList[9+i];
		rnaDihedImpD4D5Shift[i] = pList[15+i];
		rnaDihedD2D3D4Shift[i] = pList[21+i];
	}
	for(int i=0;i<4;i++){
		rnaRiboseRotamerShift[i] = pList[27+i];
	}
	for(int i=0;i<16;i++){
		wtRiboseOxy[i] = pList[31+i];
	}

}

ForceFieldPara::~ForceFieldPara() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPmodel */
