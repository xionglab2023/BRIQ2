/*
 * BasePair6DEnergyTableCG.h
 *
 *  Created on: 2023��11��28��
 *      Author: nuc
 */

#ifndef FORCEFIELD_BASEPAIR6DENERGYTABLECG_H_
#define FORCEFIELD_BASEPAIR6DENERGYTABLECG_H_

#include "dataio/datapaths.h"
#include "forcefield/ForceFieldPara.h"
#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include "forcefield/BasePair6DEnergyTable.h"
#include <time.h>
#include <math.h>
#include <map>
#include <fstream>
#include <sstream>

namespace NSPforcefield {

using namespace std;
using namespace NSPgeometry;


class BasePair6DEnergyTableCG {
public:

	map<int, double> nbKeysEnergy[36000]; //16*2250, distacne 50 bins, dihedral 45 bins, sphere 2000*2000
	map<int, double> nnbKeysEnergy[36000];

	CsMoveTo6DKey cm2Key;
	map<int,double>::iterator it;

	double wtNb, wtNnb;

	BasePair6DEnergyTableCG(ForceFieldPara* para, bool withBinary=true, int binaryMode=2);

	double getEnergy(const LocalFrame csA, const LocalFrame csB, int typeA, int typeB, int sep, double minDistance){
		if(sep > 1 && minDistance >= 5.0) return 0;
		if(minDistance < 1.5) return 0;

		double len = csA.origin_.distance(csB.origin_);
		if(len >= 15.0)
			return 0;

	    pair<int,int> p = cm2Key.toIndexPair(csA, csB, len);
		int pairType = typeA*4+typeB;
		int mapIndex = pairType*2250 + p.first;

		if(sep == 1){
			it = nbKeysEnergy[mapIndex].find(p.second);
			if(it != nbKeysEnergy[mapIndex].end()){
				//printf("indexA: %8d indexB: %8d\n", mapIndex, p.second);
				return wtNb*it->second;
			}
			else
				return 0.0;
		}
		else {
			it = nnbKeysEnergy[mapIndex].find(p.second);
			if(it != nnbKeysEnergy[mapIndex].end()){
				//printf("indexA: %8d indexB: %8d\n", mapIndex, p.second);
				return wtNnb*it->second;
			}
			else
				return 0;
		}

	}
	virtual ~BasePair6DEnergyTableCG();
};

} /* namespace NSPforcefield */

#endif /* FORCEFIELD_BASEPAIR6DENERGYTABLECG_H_ */
