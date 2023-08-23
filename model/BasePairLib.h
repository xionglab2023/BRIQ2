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

	BaseDistanceMatrix nbDMClusterCenters[16][100]; //max cluster num = 100
	BaseDistanceMatrix nnbDMClusterCenters[16][100];

	double nbEnegy[16][100];
	double nbProportion[16][100];

	double nnbEnegy[16][100];
	double nnbProportion[16][100];

	BasePairLib();

	int getPairType(BaseDistanceMatrix dm, int typeA, int typeB, int sep); //sep: sequence separation
	double getPairEnergy(RNABase* baseA, RNABase* baseB); //base-base energy + base-O2' hbond + O2'-O2' hbond


	virtual ~BasePairLib();
};

} /* namespace NSPmodel */

#endif /* MODEL_BASEPAIRLIB_H_ */
