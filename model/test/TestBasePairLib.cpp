/*
 * TestBasePairLib.cpp
 *
 *  Created on: 2023Äê8ÔÂ23ÈÕ
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

	out.close();
	delete bpLib;
	delete atLib;


}

