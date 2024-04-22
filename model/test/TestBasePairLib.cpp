/*
 * TestBasePairLib.cpp
 *
 *  Created on: 2023��8��23��
 *      Author: nuc
 */


#include "geometry/localframe.h"
#include "geometry/TransMatrix.h"
#include "model/BaseDistanceMatrix.h"
#include "model/BasePairComposition.h"
#include "model/RiboseRotamerLib.h"
#include "model/AtomLib.h"
#include "predNA/BRNode.h"
#include <time.h>
#include <stdio.h>
#include "predNA/BRFoldingTree.h"
#include "model/StructureModel.h"
#include "forcefield/AtomicClashEnergy.h"
#include "forcefield/HbondEnergy.h"
#include "model/BasePairLib.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredna;

int main(int argc, char** argv){

	BasePairLib* bpLib = new BasePairLib();
	string path2 = "bpDensityNnb";
	BasePairLib* bpLib2 = new BasePairLib(path2);
	AtomLib* atLib = new AtomLib();
	string augc = "AUGC";
	char xx[200];
	for(int i=0;i<16;i++){
		int n = bpLib2->nnbBasePairNum[i];
		for(int j=0;j<n;j++){
			BaseDistanceMatrix dm2 = bpLib2->nnbDMClusterCenters[i][j];
			double ene = bpLib2->nnbEnergy[i][j];
			double p = bpLib2->nnbProportion[i][j];

			if(ene > -9.0) continue;

			int clusterID = bpLib->getPairType(dm2, i/4, i%4, 2);
			BaseDistanceMatrix dm1 = bpLib->nnbDMClusterCenters[i][clusterID];

			double dist = dm1.distanceTo(dm2);
			if(dist < 0.5) continue;
			sprintf(xx, "/public/home/pengx/briqx/basePair/finalBasePair/pdb/nnb/%c%c%d.pdb", augc[i/4], augc[i%4], j);
			RNAPDB pdb(string(xx), "xxxx");
			vector<RNABase*> baseList = pdb.getBaseList();
			RNABase* baseA = baseList[0];
			RNABase* baseB = baseList[1];
			if(baseA->isStackingTo(baseB, atLib)) continue;
			printf("%c%c cluster: %2d %4d ene: %7.3f p: %6.4f dist: %5.3f ", augc[i/4], augc[i%4], j, clusterID, ene, p, dm1.distanceTo(dm2));
			for(int k=0;k<bpLib->nnbBasePairNum[i];k++){
				BaseDistanceMatrix dm3 = bpLib->nnbDMClusterCenters[i][k];
				if(dm3.distanceTo(dm2) < 1.2) {
					cout << " " << k;
				}
			}
			cout << endl;
		}
	}

	/*

	AtomLib* atLib = new AtomLib();

	string listFile = "/public/home/pengx/pdbLib/rna/list_R2.7";
	string outFile = "/public/home/pengx/briqx/basePair/basePairEnergyWithOxyHbond/bp.ene2";

	ofstream out;
	out.open(outFile.c_str(), ios::out);

	ifstream file;
	file.open(listFile, ios::in);
	string line;
	char xx[200];
	int i,j, sep;
	string nbnnb;

	while(getline(file, line)){
		cout << line << endl;
		string pdbID = line.substr(0, 4);
		string pdbFile = "/public/home/pengx/pdbLib/rna/pdb/"+line;
		RNAPDB pdb = RNAPDB(pdbFile, "xxxx");
		vector<RNABase*> baseList = pdb.getValidBaseList(atLib);

		for(i=0;i<baseList.size();i++){
			RNABase* baseA = baseList[i];
			int typeA = baseA->baseTypeInt;
			if(typeA < 0 || typeA > 3) continue;
			LocalFrame csA = baseA->getCoordSystem();

			for(j=i+1;j<baseList.size();j++){
				RNABase* baseB = baseList[j];
				int typeB = baseB->baseTypeInt;
				if(typeB < 0 || typeB > 3) continue;

				LocalFrame csB = baseB->getCoordSystem();

				if(j> i+1 && !baseA->contactTo(baseB))
					continue;

				if(j==i+1 && baseA->connectToNeighbor(baseB)) {
					nbnnb = "nb";
					sep = 1;
					BaseDistanceMatrix dm(csA, csB);
					double ene = bpLib->getPairEnergy(baseA, baseB);
					int type = bpLib->getPairType(dm, typeA, typeB, sep);
					double dist = bpLib->distanceToClusterCenter(dm, typeA, typeB, sep);
					if(type < 0) continue;
					if(dist > 1.2) continue;
					sprintf(xx, "%d %d %d %3d %6.3f %5.3f", typeA, typeB, sep, type, ene, dist);
					out << string(xx) << endl;
				}
				else {
					nbnnb = "nnb";
					sep = 2;
					BaseDistanceMatrix dm(csA, csB);
					double ene = bpLib->getPairEnergy(baseA, baseB);
					int type = bpLib->getPairType(dm, typeA, typeB, sep);
					double dist = bpLib->distanceToClusterCenter(dm, typeA, typeB, sep);
					if(type >=0 && dist < 1.2) {
						sprintf(xx, "%d %d %d %3d %6.3f %5.3f", typeA, typeB, sep, type, ene, dist);
						out << string(xx) << endl;
					}

					BaseDistanceMatrix dm2(csB, csA);
					double ene2 = bpLib->getPairEnergy(baseB, baseA);
					int type2 = bpLib->getPairType(dm2, typeB, typeA, sep);
					double dist2 = bpLib->distanceToClusterCenter(dm2, typeB, typeA, sep);
					if(type2 >=0 && dist2 < 1.2) {
						sprintf(xx, "%d %d %d %3d %6.3f %5.3f", typeB, typeA, sep, type2, ene2, dist2);
						out << string(xx) << endl;
					}
				}

			}
		}
	}
	delete atLib;
	out.close();
	*/


	
	delete bpLib;
	delete bpLib2;


}

