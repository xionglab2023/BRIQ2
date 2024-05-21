/*
 * TestGraph.cpp
 *
 *  Created on: 2023��12��1��
 *      Author: nuc
 */



#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>

#include "model/StructureModel.h"
#include "predNA/NuGraph.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredNA;
using namespace std;


int main(int argc, char** argv){


	clock_t start = clock();

	string inputFile = string(argv[1]);
	string output = string(argv[2]);
	int randseed = atoi(argv[3]);

	srand(randseed);

	string libType = "stat";

	BasePairLib* pairLib = new BasePairLib(libType);
	RotamerLib* rotLib = new RotamerLib();
	AtomLib* atLib = new AtomLib();

	cout << "load moveLib" << endl;
	NuPairMoveSetLibrary* moveLib = new NuPairMoveSetLibrary(libType, true, 1);
	moveLib->load();

	EdgeInformationLib* eiLib = new EdgeInformationLib();

	cout << "load energy table" << endl;
	ForceFieldPara* para = new ForceFieldPara();
	para->libType = libType;

	RnaEnergyTable* et = new RnaEnergyTable(para);
	
	et->para->kStepNum1 = 400;
	et->para->kStepNum2 = 100;
	et->para->kStepNum3 = 100;

	et->loadAtomicEnergy();

	cout << "init graph" << endl;
	NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, eiLib, et);
	graph->initForMC(inputFile);
	
	//graph->initInfo->printPDB(output);

	graph->initRandWeight();
	cout << "all edge" << endl;
	graph->printAllEdge();
	NuTree* tree = new NuTree(graph);
	graph->MST_kruskal(tree);
	cout << "tree: " << endl;
	tree->printEdges();
	cout << "adde node info" << endl;
	tree->updateNodeInfo(1.0, 1.0);
	
	cout << "edgeInfo: " << endl;
	tree->updateEdgeInfo(1.0, 1.0);
	
	for(int i=0;i<tree->geList.size();i++){
		cout << "edge: " << tree->geList[i]->indexA << "-" << tree->geList[i]->indexB << endl;
		tree->geList[i]->printPartition();
	}

	cout << "update sampling info" << endl;
	tree->updateSamplingInfo();
	//tree->printNodeInfo();


	tree->printEdgeInfo();
	
	cout << "run MC" << endl;
	graphInfo* gi = tree->runAtomicMC();
	gi->printPDB(output);
	delete gi;

	delete pairLib;
	delete rotLib;
	delete atLib;
	delete moveLib;
	delete et;
	delete tree;
	delete graph;


	clock_t end1 = clock();
	cout << "time1: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;

}


