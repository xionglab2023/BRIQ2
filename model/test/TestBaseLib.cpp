/*
 * TestBaseLib.cpp
 *
 */


#include "model/RNABaseLib.h"
#include "geometry/RMSD.h"
#include "stdlib.h"
#include <iostream>
#include "model/StructureModel.h"
#include "model/BasePairLib.h"
#include <model/BaseDistanceMatrix.h>

using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;

int main(int argc, char** argv){

	if(argc < 2){
		cout << "testBaseLib $PDB" << endl;
		exit(0);
	}
	string pdbFile = string(argv[1]);
	RNAPDB* pdb = new RNAPDB(pdbFile, "xxxx");
	vector<RNABase*> baseList = pdb->getBaseList();

	BasePairLib* bpLib = new BasePairLib();

	for(int i=0;i<baseList.size()-1;i++){
		RNABase* baseA = baseList[i];
		RNABase* baseB = baseList[i+1];
		if(!baseA->connectToNeighbor(baseB)) continue;

		LocalFrame csA = baseA->getCoordSystem();
		LocalFrame csB = baseB->getCoordSystem();
		int typeA = baseA->baseTypeInt;
		int typeB = baseB->baseTypeInt;
		BaseDistanceMatrix dm(csA, csB);

		BaseDistanceMatrix revDM(csB, csA);
		BaseDistanceMatrix revDM2 = dm.reverse();


		string augc = "AUGC";
		
		
		int pairID = bpLib->getPairType(dm, typeA, typeB, 1);
		int pairIDRev = bpLib->getPairType(revDM, typeB, typeA, -1);

		double ene1 = bpLib->getEnergy(dm, typeA, typeB, 1);
		double ene2 = bpLib->getEnergy(revDM, typeB, typeA, -1);

		double p1 = bpLib->getPairClusterProportion(dm, typeA, typeB, 1);
		double p2 = bpLib->getPairClusterProportion(revDM, typeB, typeA, -1);
		
		vector<int> neighborClusters;
		vector<double> distanceToClusterCenters;

		vector<int> neighborClusters2;
		vector<double> distanceToClusterCenters2;

		bpLib->getNeighborClusters(dm, typeA, typeB, 1, neighborClusters, distanceToClusterCenters);
		bpLib->getNeighborClusters(revDM, typeB, typeA, -1, neighborClusters2, distanceToClusterCenters2);


		double d1 = dm.distanceTo(bpLib->nbDMClusterCenters[typeA*4+typeB][pairID]);
		double d2 = revDM.distanceTo(bpLib->revNbDMClusterCenters[typeB*4+typeA][pairIDRev]);


		printf("type: %c%c seqID: %3s pairID: %4d %4d energy: %7.3f %7.3f probability: %6.4f %6.4f distance: %5.3f %5.3f neighbor: %d %d %5.3f %5.3f\n", augc[typeA], augc[typeB], baseA->baseID.c_str(), pairID, pairIDRev, ene1, ene2, p1, p2, d1, d2, neighborClusters.size(), neighborClusters2.size(), distanceToClusterCenters[0], distanceToClusterCenters2[0]);

	}

	delete bpLib;
	delete pdb;
}



