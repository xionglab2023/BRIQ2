/*
 * BasePairLib.h
 *
 *  Created on: 2023Äê8ÔÂ7ÈÕ
 *      Author: pengx
 */

#ifndef MODEL_BASEPAIRLIB_H_
#define MODEL_BASEPAIRLIB_H_

#include "geometry/xyz.h"
#include "geometry/Angles.h"
#include <vector>
#include "geometry/localframe.h"
#include "geometry/CsMove.h"
#include "model/AtomLib.h"
#include "model/ResName.h"
#include "model/StructureModel.h"
#include "model/BaseDistanceMatrix.h"
#include "tools/StringTool.h"
#include "dataio/datapaths.h"


namespace NSPmodel {


class BasePairLib {
public:

	int nbBasePairNum[16];
	int nnbBasePairNum[16];

	BaseDistanceMatrix nbDMClusterCenters[16][200]; //max cluster num < 200
	BaseDistanceMatrix nnbDMClusterCenters[16][200];

	double nbEnegy[16][200];
	double nbEnergyWithOxy[16][200];
	double nbProportion[16][200];

	double nnbEnegy[16][200];
	double nnbEnergyWithOxy[16][200];
	double nnbProportion[16][200];


	BasePairLib();

	int getPairType(BaseDistanceMatrix dm, int typeA, int typeB, int sep); //sep: sequence separation
	double distanceToClusterCenter(BaseDistanceMatrix dm, int typeA, int typeB, int sep);
	double getEnergy(BaseDistanceMatrix dm, int typeA, int typeB, int sep);
	double getPairEnergy(RNABase* baseA, RNABase* baseB); //base-base energy + base-O2' hbond + O2'-O2' hbond

	double getEnergyWithOxy(BaseDistanceMatrix dm, int typeA, int typeB, int sep);
	double getEnergyWithOxy(int clusterID, int typeA, int typeB, int sep){
		if(sep == 1)
			return nbEnergyWithOxy[typeA*4+typeB][clusterID];
		else if(sep == 2)
			return nnbEnergyWithOxy[typeA*4+typeB][clusterID];
		else if(sep == -1)
			return nbEnergyWithOxy[typeB*4+typeA][clusterID];
		else
			return 0.0;
	}

	virtual ~BasePairLib();
};

} /* namespace NSPmodel */

#endif /* MODEL_BASEPAIRLIB_H_ */
