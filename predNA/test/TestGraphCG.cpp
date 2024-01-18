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



	string inputFile = string(argv[1]);
	string output = string(argv[2]);
	int randseed = atoi(argv[3]);

	srand(randseed);

	BasePairLib* pairLib = new BasePairLib();
	RotamerLib* rotLib = new RotamerLib();
	AtomLib* atLib = new AtomLib();
	NuPairMoveSetLibrary* moveLib = new NuPairMoveSetLibrary();
	RnaEnergyTable* et = new RnaEnergyTable();
    et->loadCoarseGrainedEnergy();

	cout << "init graph" << endl;
	NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, et);
	graph->initRandWeight();
	cout << "all edge" << endl;
	graph->printAllEdge();
	NuTree* tree = new NuTree(graph);
	graph->MST_kruskal(tree);
	cout << "tree: " << endl;
	tree->printEdges();
	cout << "adde node info" << endl;
	tree->updateNodeInfoCG();
    
	for(int i=0;i<graph->seqLen;i++){
		graph->allNodes[i]->printNodeInfo();
	}

	cout << "edgeInfo: " << endl;
	tree->updateEdgeInfoCG();

	for(int i=0;i<tree->geList.size();i++){
		cout << "edge: " << tree->geList[i]->indexA << "-" << tree->geList[i]->indexB << endl;
		tree->geList[i]->printPartition();
	}


	cout << "update sampling info" << endl;
	tree->updateSamplingInfo();
	tree->printNodeInfo();

	cout << "run MC" << endl;

	clock_t start = clock();

	tree->runCoarseGrainedMC(output);
	clock_t end1 = clock();
	cout << "time1: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;

	delete pairLib;
	delete rotLib;
	delete atLib;
	delete moveLib;
	delete et;
	delete tree;
	delete graph;



}


