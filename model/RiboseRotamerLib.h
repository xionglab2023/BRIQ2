/*
 * RiboseRotamerLib.h
 *
 *  Created on: 2022��9��6��
 *      Author: pengx
 */

#ifndef MODEL_RIBOSEROTAMERLIB_H_
#define MODEL_RIBOSEROTAMERLIB_H_

#include <time.h>
#include "model/StructureModel.h"
#include "model/RiboseRotamer.h"
#include "forcefield/XPara.h"

namespace NSPmodel {

using namespace NSPforcefield;

class RiboseRotamerLib {
public:
	RiboseRotamer* rotLib[8][1500];

	RiboseRotamerLib();
	RiboseRotamerLib(XPara* para);

	RiboseRotamer* getLowestEnergyRotamer(int type){
		if(type == 0)
			return rotLib[0][264];
		else if(type == 1)
			return rotLib[1][361];
		else if(type == 2)
			return rotLib[2][78];
		else if(type == 3)
			return rotLib[3][398];
		else if(type == 4)
			return rotLib[4][1375];
		else if(type == 5)
			return rotLib[5][1000];
		else if(type == 6)
			return rotLib[6][671];
		else if(type == 7)
			return rotLib[7][567];
		else
			return NULL;
	}

	RiboseRotamer* getRandomRotamerLv1(int baseType){
		return rotLib[baseType][rand()%1500];
	}

	RiboseRotamer* getFlipRotamer(int baseType){
		/*
		 * RNA
		 * type0: imp < 0 && chi > -55 && chi < 160    600
		 * type1: imp < 0 && (chi < -55 || chi > 160)  300
		 * type2: imp > 0 && chi > -55 && chi < 160    300
		 * type3: imp > 0 && (chi < -55 || chi > 160)  300
		 * DNA
		 * type4: imp < 0 && chi > -55 && chi < 160    300
		 * type5: imp < 0 && (chi < -55 || chi > 160)   75
		 * type6: imp > 0 && chi > -55 && chi < 160   1050
		 * type7: imp > 0 && (chi < -55 || chi > 160)   75
		 */

		int r1 = rand()%10;
		if(r1 < 2)
			return rotLib[baseType][rand()%600];
		else if(r1 < 6)
			return rotLib[baseType][rand()%300 + 600];
		else if(r1 < 7)
			return rotLib[baseType][rand()%300 + 900];
		else
			return rotLib[baseType][rand()%300 + 1200];
	}

	RiboseRotamer* getNearestRotamer(RNABase* base){

		int type = base->baseTypeInt;
		double minDist = 9999.9;
		RiboseRotamer* rot = new RiboseRotamer(base);
		RiboseRotamer* nearest;
		for(int i=0;i<1500;i++){
			double d = rotLib[type][i]->distanceTo(rot);
			if(d < minDist){
				minDist = d;
				nearest = rotLib[type][i];
			}
		}
		delete rot;
		return nearest;
	}

	virtual ~RiboseRotamerLib();
};

} /* namespace NSPmodel */

#endif /* MODEL_RIBOSEROTAMERLIB_H_ */
