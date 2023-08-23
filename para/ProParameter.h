/*
 * ProParameter.h
 *
 */

#ifndef PARA_PROPARAMETER_H_
#define PARA_PROPARAMETER_H_

#include "model/AtomLib.h"
#include "tools/StringTool.h"

namespace NSPpara {

using namespace NSPmodel;
using namespace NSPtools;

class ProParameter {
public:
	AtomLib* atLib;
	string curve;

	double vdwRange;
	double wCore;
	double wSurf;
	double saiD0;
	double saiSigma;
	double wS1;
	double wS2;
	double wS2Sep[5];
	double wS245[45];
	double wtRot;
	double dsPara[4];
	double vdwSepWeight[5];
	double paraVdw[525]; //shiftS, lamdaS, shiftL, lamdaL, wellDepth: 105*5
	double vdwNbRescale[5];
	double ref[20];

	double polarSepWeight[5];
	double paraPolar[238][6]; //wdS, denRangeS, kPosS,  wdL,  denRangeL, kPosL
	double polarDamping[28];

	double vdwShiftS[167][167];
	double vdwLamdaS[167][167];
	double vdwShiftL[167][167];
	double vdwLamdaL[167][167];
	double vdwWd[167][167];
	double vdwNbD0[167][167];

	ProParameter(); //abacus2 energy
	//ProParameter(int type, const string& file, const string& filePolar);
	ProParameter(bool designTag, const string& file, const string& filePolar);
	ProParameter(const string& file);
	ProParameter(const string& fileVdw, const string& filePolar);
	void updateVdwWellDepth();
	void updateVdwWellDepthDesign();

	inline double energyRescale(double e){
		if(e < 2)
			return e;
		else
			return 2*sqrt(e*2) - 2;
	}
	virtual ~ProParameter();
};

} /* namespace NSPmodel */

#endif /* PARA_PROPARAMETER_H_ */
