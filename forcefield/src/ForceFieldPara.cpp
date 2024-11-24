/*
 * ForceFieldPara.cpp
 *
 *  Created on: 2023��4��25��
 *      Author: pengx
 */

#include "forcefield/ForceFieldPara.h"

namespace NSPforcefield {

ForceFieldPara::ForceFieldPara() {
	// TODO Auto-generated constructor stub
	this->clashLamdaNb = 2.8;
	this->clashLamdaNnb = 2.8;

	this->hbLamda1 = 0.15;
	this->hbLamda2 = 0.3;
	this->kHbOri = 1.0;
	this->wtHb = 1.8;

	this->rnaKAng = 0.1;
	this->rnaKBond = 9.0;

	this->wtPho = 0.7;
	this->wtRibose = 1.0;

	for(int i=0;i<6;i++){
		rnaDihedImpD1D2Shift[i] = 0.0;
		rnaDihedImpD4D5Shift[i] = 0.0;
		rnaDihedD2D3D4Shift[i] = 0.0;
	}

	rnaDihedImpD4D5Shift[0] = -1.0;
	rnaDihedImpD4D5Shift[3] = -1.0;

	for(int i=0;i<16;i++){
		rnaRiboseRotamerShift[i] = 0.0;
	}


	rnaRiboseRotamerShift[0] = 1.3;   rnaRiboseRotamerShift[1] = -1.8;
	rnaRiboseRotamerShift[4] = 1.4;   rnaRiboseRotamerShift[5] = -2.0;
	rnaRiboseRotamerShift[8] = 2.3;   rnaRiboseRotamerShift[9] = -1.5;
	rnaRiboseRotamerShift[12] = 2.4;  rnaRiboseRotamerShift[13] = -1.5;


	for(int i=0;i<16;i++){
		wtRiboseOxy[i] = 1.5;
	}

	this->bwTag = "bw5";
	this->libType = "stat"; //"stat" or "xtb" or "adj"
	this->cgEneType = "stat"; //"stat" or "xtb"

	this->wtBp1 = 2.0;
	this->wtBp2 = 2.0;

	this->wtO4O2C2Nb = 0.7;
	this->wtO4O2C2Nnb = 0.3;

	this->wtClash = 1.0;

	this->T0 = 15.0;
	this->T1 = 1.5;
	this->T2 = 0.15;
	this->T3 = 0.015;

	this->kStepNum1 = 500;
	this->kStepNum2 = 200;
	this->kStepNum3 = 100;

	this->kStepNum1CG = 200;
	this->kStepNum2CG = 80;
	this->kStepNum3CG = 50;

	this->kNodeFreq = 1.0;

	for(int i=0;i<16;i++){
		for(int j=0;j<6;j++){
			this->nbPairEnergyRescale[i][j] = 1.0;
		}
	}

	this->clashRescale = 0.1;
	this->connectRescale = 0.4;
	this->rescaleIncreaseFactor = 1.05;

	this->phoRep = 2.0;

	addSubEnergyPairtype = -1;
	addSubEnergyClusterID = -1;
	addSubEnergyWeight = 0.0;

	this->withRandomInit = true;

}

ForceFieldPara::ForceFieldPara(const string& paraFile){
	this->clashLamdaNb = 2.8;
	this->clashLamdaNnb = 2.8;

	this->hbLamda1 = 0.15;
	this->hbLamda2 = 0.3;
	this->kHbOri = 1.0;
	this->wtHb = 1.8;

	this->rnaKAng = 0.1;
	this->rnaKBond = 9.0;

	this->wtPho = 0.7;
	this->wtRibose = 1.0;

	for(int i=0;i<6;i++){
		rnaDihedImpD1D2Shift[i] = 0.0;
		rnaDihedImpD4D5Shift[i] = 0.0;
		rnaDihedD2D3D4Shift[i] = 0.0;
	}

	rnaDihedImpD4D5Shift[0] = -1.0;
	rnaDihedImpD4D5Shift[3] = -1.0;

	for(int i=0;i<16;i++){
		rnaRiboseRotamerShift[i] = 0.0;
	}


	rnaRiboseRotamerShift[0] = 1.3;   rnaRiboseRotamerShift[1] = -1.8;
	rnaRiboseRotamerShift[4] = 1.4;   rnaRiboseRotamerShift[5] = -2.0;
	rnaRiboseRotamerShift[8] = 2.3;   rnaRiboseRotamerShift[9] = -1.5;
	rnaRiboseRotamerShift[12] = 2.4;  rnaRiboseRotamerShift[13] = -1.5;


	for(int i=0;i<16;i++){
		wtRiboseOxy[i] = 1.5;
	}

	this->bwTag = "bw5";
	this->libType = "stat"; //"stat" or "xtb" or "adj"
	this->cgEneType = "stat"; //"stat" or "xtb"

	this->wtBp1 = 2.0;
	this->wtBp2 = 2.0;

	this->wtO4O2C2Nb = 0.7;
	this->wtO4O2C2Nnb = 0.3;

	this->wtClash = 1.0;

	this->T0 = 15.0;
	this->T1 = 1.5;
	this->T2 = 0.15;
	this->T3 = 0.015;

	this->kStepNum1 = 1280;
	this->kStepNum2 = 1280;
	this->kStepNum3 = 1280;

	this->kStepNum1CG = 200;
	this->kStepNum2CG = 80;
	this->kStepNum3CG = 50;

	this->kNodeFreq = 1.0;

	for(int i=0;i<16;i++){
		for(int j=0;j<6;j++){
			this->nbPairEnergyRescale[i][j] = 1.0;
		}
	}

	this->clashRescale = 0.1;
	this->connectRescale = 0.4;
	this->rescaleIncreaseFactor = 1.05;

	this->phoRep = 2.0;

	addSubEnergyPairtype = -1;
	addSubEnergyClusterID = -1;
	addSubEnergyWeight = 0.0;

	this->withRandomInit = true;

}

ForceFieldPara::~ForceFieldPara() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPmodel */
