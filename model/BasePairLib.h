/*
 * BasePairLib.h
 *
 *  Created on: 2023��8��7��
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
	int nbContactBasePairNum[16];
	int nbNonContactBasePairNum;
	int nnbBasePairNum[16];

	BaseDistanceMatrix nbDMClusterCenters[16][3000]; //max cluster num < 3000
	BaseDistanceMatrix revNbDMClusterCenters[16][3000]; 
	BaseDistanceMatrix nnbDMClusterCenters[16][3000];

	/*
	 * all neighbor pairs
	 */
	double nbEnergy[16][3000];
	double nbEnergyWithOxy[16][3000];
	double nbProportion[16][3000];

	int reversePairClusterID[16][3000];

	/*
	 * non-neighbor pairs
	 */
	double nnbEnergy[16][3000];
	double nnbEnergyWithOxy[16][3000];
	double nnbProportion[16][3000];

	string libType;

	BasePairLib(const string& libType = "stat");

	int getPairType(BaseDistanceMatrix& dm, int typeA, int typeB, int sep, double ddmCutoff=1.2); //sep: sequence separation

	int getNeighborPairFirstFiveClusterID(BaseDistanceMatrix& dm, int typeA, int typeB);
	void getNeighborClusters(BaseDistanceMatrix& dm, int typeA, int typeB, int sep, vector<int>& neighborClusters, vector<double>& distanceToClusterCenters, double distanceCutoff = 1.2);
	void getNearestCluster(BaseDistanceMatrix& dm, int typeA, int typeB, int sep, vector<int>& neighborClusters, vector<double>& distanceToClusterCenters);
	double getPairClusterProportion(BaseDistanceMatrix& dm, int typeA, int typeB, int sep);
	double distanceToClusterCenter(BaseDistanceMatrix& dm, int typeA, int typeB, int sep);
	double getEnergy(BaseDistanceMatrix& dm, int typeA, int typeB, int sep);
	
	double getEnergy(int clusterID, int typeA, int typeB, int sep){
		if(clusterID < 0)
			return 0.0;
		if(sep == 1)
			return nbEnergy[typeA*4+typeB][clusterID];
		else if(sep == 2)
			return nnbEnergy[typeA*4+typeB][clusterID];
		else if(sep == -1)
			return nbEnergy[typeB*4+typeA][clusterID];
		else
			return 0.0;
	}

	double getDMDistance(int clusterID1, int clusterID2, int typeA, int typeB, int sep){
		if(clusterID1 < 0 && clusterID2 < 0)
			return 0.0;
		else if(clusterID1 < 0 || clusterID2 < 0)
			return 3.0;

		if(sep == 1)
			return nbDMClusterCenters[typeA*4+typeB][clusterID1].distanceTo(nbDMClusterCenters[typeA*4+typeB][clusterID2]);
		else if(sep == 2)
			return nnbDMClusterCenters[typeA*4+typeB][clusterID1].distanceTo(nnbDMClusterCenters[typeA*4+typeB][clusterID2]);
		else if(sep == -1)
			return nbDMClusterCenters[typeB*4+typeA][clusterID1].distanceTo(nbDMClusterCenters[typeB*4+typeA][clusterID2]);
		else
			return 0.0;
	}

	double getEnergyWithOxy(int clusterID, int typeA, int typeB, int sep){
		if(clusterID < 0)
			return 0.0;
		if(sep == 1)
			return nbEnergyWithOxy[typeA*4+typeB][clusterID];
		else if(sep == 2)
			return nnbEnergyWithOxy[typeA*4+typeB][clusterID];
		else if(sep == -1)
			return nbEnergyWithOxy[typeB*4+typeA][clusterID];
		else
			return 0.0;
	}

	double getPairEnergy(RNABase* baseA, RNABase* baseB); //base-base energy + base-O2' hbond + O2'-O2' hbond

	virtual ~BasePairLib();
};

} /* namespace NSPmodel */

#endif /* MODEL_BASEPAIRLIB_H_ */
