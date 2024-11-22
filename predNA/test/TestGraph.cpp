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

	EdgeMoveClustersLib* eiLib = new EdgeMoveClustersLib();

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

	cout << "int for MC" << endl;
	graph->initForMC(inputFile);
	
	cout << "graph generate random partition" << endl;
	graph->generateRandomEdgePartition(2);
	SamplingGraph* sg = new SamplingGraph(graph);

	cout << "SG update partition" << endl;
    sg->updatePartitionInfo();

	cout << "SG update sampling" << endl;
    sg->updateSamplingInfo();
	
	cout << "final partition: " << sg->geList.size() << endl;
	for(int i=0;i<sg->geList.size();i++){
		cout << "p edge: " << sg->geList[i]->indexA << "-" << sg->geList[i]->indexB << " " << sg->geList[i]->epList.size() << endl;
	}

	sg->printEdgeInfo();
	cout << endl;
	sg->printEdges();

	cout << "before MC" << endl;
	sg->checkCsMove();
	
	cout << "run cg mc" << endl;
	sg->runCoarseGrainedMC(output);

	graph->printEdgeClusterRegister();

	//graphInfo* gi = sg->runAtomicMC();
	//gi->printPDB(output);
	//delete gi;

	delete pairLib;
	delete rotLib;
	delete atLib;
	delete moveLib;
	delete et;
	delete sg;
	delete graph;


	clock_t end1 = clock();
	cout << "time1: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;

}


