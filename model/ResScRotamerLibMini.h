/*
 * ResScRotamerLibMini.h
 *
 *  Created on: 2022Äê7ÔÂ16ÈÕ
 *      Author: pengx
 */

#ifndef MODEL_RESSCROTAMERLIBMINI_H_
#define MODEL_RESSCROTAMERLIBMINI_H_
#include "geometry/CsMove.h"
#include "geometry/TransMatrix.h"
#include "model/ResScRotamer.h"
#include "dataio/datapaths.h"
#include "tools/StringTool.h"
#include <vector>
#include <fstream>
#include <time.h>

namespace NSPmodel {

class ResScRotamerLibMini {
public:

	AtomLib* atLib;
	vector<ResScRotamer*> rotList[20];
	vector<int> rotNum;

	//backbone rotamer index: amino acid specific

	double eRot[20][1000][80];

	ResScRotamerLibMini();

	int getRotamerID(ResScRotamer* rot);

	double getEnergy(int bbIndex, ResScRotamer* rot){

		if(bbIndex < 0 || bbIndex > 999) {
			cout << "invalid bbIndex" << endl;
			exit(0);
		}
		int type = rot->aaType;
		if(type < 0 || type >= 20) {
			cout << "invalid aa type " << type << endl;
		}

		return eRot[type][bbIndex][rot->rotID];
	}

	virtual ~ResScRotamerLibMini();
};

}
#endif /* MODEL_RESSCROTAMERLIBMINI_H_ */
