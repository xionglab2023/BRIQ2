/*
 * XPara.h
 *
 *  Created on: 2023Äê5ÔÂ6ÈÕ
 *      Author: nuc
 */

#ifndef FORCEFIELD_XPARA_H_
#define FORCEFIELD_XPARA_H_

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "tools/StringTool.h"

namespace NSPforcefield {

using namespace std;

class XPara {

public:
	double wtBp1;
	double wtBp2;
	double wtBp3;

	double wtRotamer;
	//double wtBaseOxygen;
	double kBond;
	double kAng;
	double wtDihed;
	double wtPho;
	double wtConnect;

	string dihedEneType;
	string rotEneType;
	string atomicEneType;
	string clashType;

	double wtClash;
	double lamdaClash;
	double bbClash;

	double T0;
	double T1;
	double anneal;

	double initShift;
	double dShift;

	int stepNum;

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

	/*
	int smcModelNum;
	double smcT0;
	double smcT1;
	int smcStepNum;
	int smcClusterNum;

	double omcT0;
	double omcT1;
	int omcStepNum;
	*/

	XPara();
	XPara(const string& paraFile);

	virtual ~XPara();
};

inline double energyRescale(double e){
	if(e < 0.81)
		return e;
	else
		return 1.8*sqrt(e) - 0.81;
}

inline double energyRescale8(double e){
	//return 0;
	if(e < 8.0)
		return e;
	else {
		return 8.0*log(e) - 8.635532333;
	}
}

} /* namespace NSPforcefield */


#endif /* FORCEFIELD_XPARA_H_ */
