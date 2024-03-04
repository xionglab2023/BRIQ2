/*
 * ForceFieldPara.h
 *
 *  Created on: 2023��4��25��
 *      Author: pengx
 */

#ifndef FORCEFIELD_FORCEFIELDPARA_H_
#define FORCEFIELD_FORCEFIELDPARA_H_

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "tools/StringTool.h"

namespace NSPforcefield {

using namespace std;

class ForceFieldPara {
public:

	double clashLamdaNb;
	double clashLamdaNnb;

	double hbLamda1;
	double hbLamda2;
	double kHbOri;
	double wtHb;

	double rnaKBond;
	double rnaKAng;

	double wtPho;
	double wtRibose;

	double rnaDihedImpD1D2Shift[6];
	double rnaDihedImpD4D5Shift[6];
	double rnaDihedD2D3D4Shift[6];

	double rnaRiboseRotamerShift[16];

	double wtRiboseOxy[16];

	double wtBp1;
	double wtBp2;

	double wtO4O2C2Nb;
	double wtO4O2C2Nnb;

	double wtClash;
	double lamdaClash;
//	double bbClash;
	string bwTag;

	double T0;
	double T1;
	double T2;
	double T3;

	int kStepNum1;
	int kStepNum2;
	int kStepNum3;
	double kNodeFreq;


	double anneal;

	double initShift;
	double dShift;



	int outFreq;

	double initConnectWT;
	double initClashWT;
	double connectWTFactor;
	double clashWTFactor;

	bool loopRiboConnectMove;
	bool ctRandMove;
	bool f3Move;
	bool singleBaseMove;
	bool reverseRotMove;

	string spType;
	double phoRep;


	ForceFieldPara();
	ForceFieldPara(const string& paraFile);
	virtual ~ForceFieldPara();
};

} /* namespace NSPmodel */

#endif /* FORCEFIELD_FORCEFIELDPARA_H_ */
