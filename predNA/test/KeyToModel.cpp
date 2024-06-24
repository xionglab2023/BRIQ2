/*
 * KeyToModel.cpp
 *
 *  Created on: 2024,6,22
 *      Author: nuc
 */



#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include "tools/CmdArgs.h"
#include "model/StructureModel.h"
#include "predNA/NuGraph.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredNA;
using namespace std;

void printHelp(){
	cout << "key2pdb -in $INPUT -out $OUTPUTPDB -ene $OUTENEFILE -seed $RANDSEED"
}


int main(int argc, char** argv){
	clock_t start = clock();
    CmdArgs cmdArgs{argc, argv};

	if(argc == 1 || cmdArgs.specifiedOption("-h")){
		printHelp();
	} 

    string inputFile = cmdArgs.getValue("-in");
    string output = cmdArgs.getValue("-out");
    string eneout = cmdArgs.getValue("-ene");
    string seed = cmdArgs.getValue("-seed");

	srand(atoi(seed.c_str()));

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
	

	cout << "update sampling info" << endl;
	tree->updateSamplingInfo();

	
	cout << "run MC" << endl;
	graphInfo* gi = tree->runAtomicMC();
	//graph->printEnergy();

	gi->printPDB(output);

    gi->printDetailEnergy(eneout, pairLib, atLib, et);
    
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


