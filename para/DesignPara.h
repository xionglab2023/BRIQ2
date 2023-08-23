/*
 * DesignPara.h
 *
 *  Created on: 2022Äê3ÔÂ10ÈÕ
 *      Author: pengx
 */

#ifndef PARA_DESIGNPARA_H_
#define PARA_DESIGNPARA_H_

#include "model/AtomLib.h"
#include "tools/StringTool.h"

namespace NSPpara {

using namespace NSPmodel;

class DesignPara {

public:
	AtomLib* atLib;


	double vdwSepWeight[5];
	string curve;
	double vdwRange;
	double wd0;
	double vdwWdRescale[60]; //15*4
	double polarWdRescale[4];

	double wtRot;

	double shiftS;
	double shiftL;
	double lamdaS;
	double lamdaL;

	double wdHb;

	double dsPara[4];
	double wS1;
	double wS245[45];



	double vdwRescale1[167][167];
	double vdwRescale2[167][167];
	double vdwRescale3[167][167];
	double vdwRescale4[167][167];

	double vdwNbD0[167][167];

	double polarKS;
	double polarKL;
	double polarCutS;
	double polarCutL;

	double polarWd[70]; //polar group 14 * polar type 5
	double ref[60];

	DesignPara(const string& file);
	void updateVdwWellDepthDesign();


	virtual ~DesignPara();
};

} /* namespace NSPpara */

#endif /* PARA_DESIGNPARA_H_ */
