/*
 * TestGraph.cpp
 *
 *  Created on: 2024,3,4
 *      Author: pengx
 */

#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include "model/StructureModel.h"
#include "predNA/NuGraph.h"
#include "tools/CmdArgs.h"
#include "tools/ThreadPool.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredNA;
using namespace std;
using namespace NSPtools;
using namespace NSPthread;


int runCGMC(NuPairMoveSetLibrary* moveLib, RnaEnergyTable* et, const string& inputFile, const string& outFile, int randSeed){

 	srand(randSeed);
	BasePairLib* pairLib = new BasePairLib();
	RotamerLib* rotLib = new RotamerLib();
	AtomLib* atLib = new AtomLib();
    et->loadCoarseGrainedEnergy();

	NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, et);
	graph->initForCGMC(inputFile);
	graph->initRandWeight();
	graph->printAllEdge();
	NuTree* tree = new NuTree(graph);
	graph->MST_kruskal(tree);
	tree->printEdges();
	tree->updateNodeInfoCG();
    
	for(int i=0;i<graph->seqLen;i++){
		graph->allNodes[i]->printNodeInfo();
	}

	tree->updateEdgeInfoCG();

	for(int i=0;i<tree->geList.size();i++){
		tree->geList[i]->printPartition();
	}

	tree->updateSamplingInfo();
	tree->printNodeInfo();

	clock_t start = clock();

	tree->runCoarseGrainedMC(outFile);
	clock_t end1 = clock();
	cout << "time1: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;

	delete pairLib;
	delete rotLib;
	delete atLib;
	delete tree;
	delete graph;

    return 0;
}

void printHelp(){
    cout << "testGraphCG -in $INPUTFILE -out $OUTPUTFILE -mp 64" << endl;
}

int main(int argc, char** argv){

    if(argc == 1) {
        printHelp();
        return EXIT_SUCCESS;
    }
    CmdArgs cmdArgs{argc, argv};


	NuPairMoveSetLibrary* moveLib = new NuPairMoveSetLibrary(true, 1);
	moveLib->load();

	RnaEnergyTable* et = new RnaEnergyTable();
	et->loadAtomicEnergy();

    string inputFile = cmdArgs.getValue("-in");
    string outputFile = cmdArgs.getValue("-out");

    string outEndTag = outputFile.substr(outputFile.length()-3, 3);

    if(outEndTag != "pdb") {
        cout << "output file should be end with .pdb" << endl;
        exit(0);
    }

    int mp = atoi(cmdArgs.getValue("-mp").c_str());
    int startID = 0;
    if(cmdArgs.specifiedOption("-id")) {
        startID = atoi(cmdArgs.getValue("-id").c_str());
    }

	clock_t start = clock();

    shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
    size_t jid = 0;
    char xx[200];
    for(int i=startID;i<startID+mp;i++) {
        shared_ptr<IntFuncTask> request(new IntFuncTask);
        sprintf(xx, "%s-%d.pdb", outputFile.substr(0, outputFile.length()-4).c_str(), i);
        string outFile2 = string(xx);
        request->asynBind(runCGMC, moveLib, et, inputFile, outFile2, time(0)+i);
        jid++;
        thrPool->addTask(request);
    }
    while(true) {
        sleep(1);
        if(thrPool->getTaskCount() == 0) {
        break;
        }
    }

    clock_t end1 = clock();
	cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;

    delete moveLib;
    delete et;
}


