/*
 * ResBBRotamerLib.h
 *
 */

#ifndef MODEL_RESBBROTAMERLIB_H_
#define MODEL_RESBBROTAMERLIB_H_

#include "geometry/CsMove.h"
#include "geometry/TransMatrix.h"
#include "dataio/datapaths.h"
#include "tools/StringTool.h"
#include <vector>
#include <fstream>
#include <time.h>
#include "model/ResBBRotamer.h"
#include "model/AtomLib.h"

namespace NSPmodel {


class ResBBRotamerLib {
public:

	AtomLib* atLib;

	ResBBRotamer* allRotLib1k[20][1000];

	int neighborLib1k[1000][50]; //nearest 50 rotamers in lib1k
	int neighborLib1w[1000][12]; //nearest 12 rotamers in lib1w

	ResBBRotamer* allRotLib1w[20][10000];

	int index1wTo1k[10000];

	ResBBRotamerLib();

	ResBBRotamer* getRandomRotamer(int aaType){
		return allRotLib1k[aaType][rand()%10000];
	}

	ResBBRotamer* getRandomNearbyRotamerP05(ResBBRotamer* initRot){
		//nearest 5% rotamers
		int i1 = neighborLib1k[initRot->index1K][rand()%50];
		int i2 = neighborLib1w[i1][rand()%12];

		if(initRot->index1W == i2)
			return getRandomNearbyRotamerP05(initRot);
		else
			return allRotLib1w[initRot->aaType][i2];
	}


	ResBBRotamer* getRandomNearbyRotamerP01(ResBBRotamer* initRot){
		//nearest 1% rotamers
		int i1 = neighborLib1k[initRot->index1K][rand()%10];
		int i2 = neighborLib1w[i1][rand()%12];

		if(initRot->index1W == i2)
			return getRandomNearbyRotamerP05(initRot);
		else
			return allRotLib1w[initRot->aaType][i2];
	}

	int getRotamerIndex1K(ResBBRotamer* bbRot){
		double minD = 99999.9;
		int minIndex = 0;
		double d;
		for(int i=0;i<1000;i++){
			d = bbRot->distanceTo(allRotLib1k[bbRot->aaType][i]);
			if(d < minD){
				minD = d;
				minIndex = i;
			}
		}
		return minIndex;
	}

	int getRotamerIndex1W(ResBBRotamer* bbRot) {
		double minD = 99999.9;
		int minIndex = 0;
		double d;
		for(int i=0;i<10000;i++){
			d = bbRot->distanceTo(allRotLib1w[bbRot->aaType][i]);
			if(d < minD){
				minD = d;
				minIndex = i;
			}
		}
		return minIndex;
	}

	virtual ~ResBBRotamerLib();
};
}
#endif /* MODEL_RESBBROTAMERLIB_H_ */
