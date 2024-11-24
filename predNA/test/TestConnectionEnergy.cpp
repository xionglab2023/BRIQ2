/*
 * TestConnectionEnergy.cpp
 *
 */

#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>

#include "model/RotamerLib.h"
#include "model/StructureModel.h"
#include "predNA/BRNode.h"
#include "predNA/BRFoldingTree.h"
#include "predNA/MCRun.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredNA;
using namespace std;


int main(int argc, char** argv){
	clock_t start = clock();
	string inputFile = string(argv[1]);

	RnaEnergyTable et;
	//RiboConnectToPO3 et;

	RNAPDB pdb(inputFile, "xxxx");
	vector<RNABase*> baseList = pdb.getBaseList();
	RotamerLib* rotLib = new RotamerLib();
	vector<BRNode*> nodes;
	RiboseRotamer* ribRot;
	for(int i=0;i<baseList.size();i++){
		BRNode* node = new BRNode(baseList[i]->baseTypeInt, i, rotLib);
		node->seqID = i;
		nodes.push_back(node);
	}

	for(int i=0;i<baseList.size()-1;i++){
		BRNode* nodeA = nodes[i];
		BRNode* nodeB = nodes[i+1];
		et.pb->buildPhosphate(nodeA->riboseConf, nodeB->riboseConf, nodeA->phoConf);
		nodeA->phoConfTmp->copyValueFrom(nodeA->phoConf);
		cout << "node: " << i << " ";
		cout << nodeA->phoConf->ene << endl;

	}

	for(int i=0;i<baseList.size()-1;i++){
		BRNode* nodeA = nodes[i];

		for(int j=i+1;j<baseList.size();j++){
			BRNode* nodeB = nodes[j];
			double e;
			if(j==i+1){
				e = getBaseBaseEnergy(nodeA, nodeB, 1, &et, false);
			}
			else if(j == i+2)
				e = getBaseBaseEnergy(nodeA, nodeB, 2, &et, false);
			else {
				e = getBaseBaseEnergy(nodeA, nodeB, 3, &et, false);
			}
			printf("%-2d %-2d %9.3f\n",i,j,e);
		}
	}


}

