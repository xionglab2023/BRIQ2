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
#include "predNA/NuSampling.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredNA;
using namespace std;
using namespace NSPtools;
using namespace NSPthread;


class cgResult{
public:
    map<string, double> keyMap;

    cgResult(){
    }

    void mergeResult(cgResult* other){
       map<string,double>::iterator it, it2;
       for(it = other->keyMap.begin();it != other->keyMap.end();++it){
            it2 = keyMap.find(it->first);
            if(it2 != keyMap.end()){
                if(it->second < it2->second)
                    keyMap[it->first] = it->second;
                else
                    keyMap[it->first] = it2->second;
            }
            else {
                keyMap[it->first] = it->second;
            }
       }
    }

    void clear(){
        this->keyMap.clear();
    }

    void print(){
        map<string,double>::iterator it, it2;
        for(it = keyMap.begin();it!=keyMap.end();++it){
            cout << it->first << " " << it->second << endl;
        }
    }
};

int runCGMC(NuPairMoveSetLibrary* moveLib, RnaEnergyTable* et, const string& inputFile, const string& outFile, cgResult* result, int randSeed){

 	srand(randSeed);
	BasePairLib* pairLibXtb = new BasePairLib("xtb");
    BasePairLib* pairLibStat = new BasePairLib("stat");
	RotamerLib* rotLib = new RotamerLib();
	AtomLib* atLib = new AtomLib();
    EdgeInformationLib* eiLib = new EdgeInformationLib();


    cout << "init graph" << endl;
	NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLibXtb, pairLibStat, moveLib, eiLib, et);

    cout << "init CGMC" << endl;
	graph->initForCGMC(inputFile);
    cout << "init rand weight" << endl;
	graph->initRandWeight();
    cout << "print edge: " << endl;
	graph->printAllEdge();
	NuTree* tree = new NuTree(graph);
	graph->MST_kruskal(tree);
	tree->printEdges();
	tree->updateNodeInfoCG(1.0, 1.0);
	tree->updateEdgeInfoCG(1.0, 1.0);
	tree->updateSamplingInfo();
	tree->printNodeInfo();
    tree->printEdgeInfo();

    string initKey = graph->toContactMapHashKeyCG();
    cout << "initKey: " << initKey << endl;
    graph->keyToContactMatrix(initKey);
	clock_t start = clock();
    cout << "run mc" << endl;
    //tree->runCoarseGrainedMC(outFile);


    NuSampling* samp = new NuSampling(graph, tree);

    samp->runCoarseGrainedMC(result->keyMap, outFile);

//    cout << "keyNum: " << result->keyMap.size() << endl;

	delete pairLibXtb;
    delete pairLibStat;
	delete rotLib;
	delete atLib;
    delete eiLib;
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


	NuPairMoveSetLibrary* moveLib = new NuPairMoveSetLibrary("stat", true, 1);
	moveLib->load();

    ForceFieldPara* para = new ForceFieldPara();
    para->libType = "stat";
    para->clashRescale = 0.1;
    para->connectRescale = 0.1;

	RnaEnergyTable* et = new RnaEnergyTable(para);
	et->loadCoarseGrainedEnergy();

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

    cgResult* resultList[mp];
    for(int i=0;i<mp;i++){
        resultList[i] = new cgResult();
    }

    for(int i=startID;i<startID+mp;i++) {
        shared_ptr<IntFuncTask> request(new IntFuncTask);
        sprintf(xx, "%s-%d.pdb", outputFile.substr(0, outputFile.length()-4).c_str(), i);
        string outFile2 = string(xx);


        request->asynBind(runCGMC, moveLib, et, inputFile, outFile2, resultList[i-startID], time(0)+i);
        jid++;
        thrPool->addTask(request);
    }

    while(true) {
        sleep(1);
        if(thrPool->getTaskCount() == 0) {
        break;
        }
    }

    for(int i=1;i<mp;i++){
        resultList[0]->mergeResult(resultList[1]);
    }

    resultList[0]->print();

    cout << "total key num: " << resultList[0]->keyMap.size() << endl;

    clock_t end1 = clock();
	cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;

    delete moveLib;
    delete et;
}


