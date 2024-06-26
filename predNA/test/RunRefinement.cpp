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

	if(argc < 4) {
		cout << "rna_refinement $input $output $mode $eneFile $randseed" << endl;
		cout << "mode: fast, normal, slow" << endl;
		exit(0);
	}

	string inputFile = string(argv[1]);
	string output = string(argv[2]);
	string mode = string(argv[3]);
	string eneFile = string(argv[4]);
	int randseed = atoi(argv[5]);

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
	

	if(mode == "fast") {
   		et->para->T0 = 0.4;
    	et->para->T1 = 0.2;
    	et->para->T2 = 0.1;
   		et->para->T3 = 0.02;
    	et->para->clashRescale = 0.5;
    	et->para->connectRescale = 0.8;
		et->para->kStepNum1 = 50;
		et->para->kStepNum2 = 30;
		et->para->kStepNum3 = 30;
		et->para->withRandomInit = false;
	}
	else if(mode == "normal") {
   		et->para->T0 = 1.0;
    	et->para->T1 = 0.2;
    	et->para->T2 = 0.05;
   		et->para->T3 = 0.01;

    	et->para->clashRescale = 0.2;
    	et->para->connectRescale = 0.5;

		et->para->kStepNum1 = 200;
		et->para->kStepNum2 = 100;
		et->para->kStepNum3 = 100;
		et->para->withRandomInit = false;
	}
	else if(mode == "slow") {
   		et->para->T0 = 5.0;
    	et->para->T1 = 0.5;
    	et->para->T2 = 0.05;
   		et->para->T3 = 0.005;

    	et->para->clashRescale = 0.1;
    	et->para->connectRescale = 0.2;

		et->para->kStepNum1 = 600;
		et->para->kStepNum2 = 300;
		et->para->kStepNum3 = 300;
		et->para->withRandomInit = true;
	}



	et->loadAtomicEnergy();

	cout << "init graph" << endl;
	NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, eiLib, et);
	graph->initForMC(inputFile);
	
	graph->initRandWeight();
	cout << "all edge" << endl;
	graph->printAllEdge();
	NuTree* tree = new NuTree(graph);
	graph->MST_kruskal(tree);
	cout << "tree: " << endl;
	tree->printEdges();
	cout << "adde node info" << endl;
	tree->updateNodeInfo(1.0, 1.0);
	for(int i=0;i<graph->seqLen;i++){
		graph->allNodes[i]->printNodeInfo();
	}

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
	gi->printDetailEnergy(eneFile, pairLib, atLib, et);
	
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


