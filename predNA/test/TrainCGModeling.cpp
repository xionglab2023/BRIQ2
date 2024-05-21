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
    map<string, double> keyToEnergy;
    map<string, double> keyToAccuracy;

    cgResult(){
    }

    void mergeResult(cgResult* other){
       map<string,double>::iterator it, it2;
       for(it = other->keyToEnergy.begin();it != other->keyToEnergy.end();++it){
            it2 = keyToEnergy.find(it->first);
            if(it2 != keyToEnergy.end()){
                if(it->second < it2->second)
                    keyToEnergy[it->first] = it->second;
                else
                    keyToEnergy[it->first] = it2->second;
            }
            else {
                keyToEnergy[it->first] = it->second;
            }
       }
    }

    void clear(){
        this->keyToEnergy.clear();
        this->keyToAccuracy.clear();
    }

    void updateAccuracy(NuGraph* natGraph){
        this->keyToAccuracy.clear();
        map<string,double>::iterator it;
        for(it = keyToEnergy.begin();it != keyToEnergy.end();++it){
            string key = it->first;
            double accuracy = natGraph->keyAccuracy(key);
            keyToAccuracy[key] = accuracy;
        }
    }

    void print(ofstream& out){
        map<string,double>::iterator it, it2;
        
        
        for(it = keyToEnergy.begin();it!=keyToEnergy.end();++it){
            if(keyToAccuracy.find(it->first) == keyToAccuracy.end())
                out << it->first << " " << it->second << " 0.0"  << endl;
            else
                out << it->first << " " << it->second << " " << keyToAccuracy[it->first] << endl;
        }
    }
};

int runCGMC(NuPairMoveSetLibrary* moveLib, RnaEnergyTable* et,EdgeInformationLib* eiLib, const string& inputFile, cgResult* result, int modelNum, int randSeed){

 	srand(randSeed);
    BasePairLib* pairLib = new BasePairLib("stat");
	RotamerLib* rotLib = new RotamerLib(et->para);
    AtomLib* atLib = new AtomLib();


    cout << "init graph" << endl;
	NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, eiLib, et);

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

    samp->runCoarseGrainedMC(result->keyToEnergy, modelNum);

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

    string inputFile = cmdArgs.getValue("-in");
    string outputFile = cmdArgs.getValue("-out");

    double rescale = atof(cmdArgs.getValue("-rescale").c_str());
    double increase = atof(cmdArgs.getValue("-increase").c_str());


	NuPairMoveSetLibrary* moveLib = new NuPairMoveSetLibrary("stat", true, 1);
	moveLib->load();

    EdgeInformationLib* eiLib = new EdgeInformationLib();
    ForceFieldPara* para = new ForceFieldPara();
    para->libType = "stat";
    para->clashRescale = rescale*0.3; //default 0.03
    para->connectRescale = rescale; //default 0.1
    para->rescaleIncreaseFactor = increase; //default 1.05


	RnaEnergyTable* et = new RnaEnergyTable(para);
	et->loadCoarseGrainedEnergy();


    ofstream out;
    out.open(outputFile.c_str(), ios::out);

    int modelNum = atoi(cmdArgs.getValue("-n").c_str());

    int mp = 2;
    int startID = 0;

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
    
        request->asynBind(runCGMC, moveLib, et, eiLib, inputFile, resultList[i-startID], modelNum, time(0)+i);
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
        resultList[0]->mergeResult(resultList[i]);
    }

    BasePairLib* pairLib = new BasePairLib("stat");
	RotamerLib* rotLib = new RotamerLib(et->para);
    AtomLib* atLib = new AtomLib();
    NuGraph* natGraph = new NuGraph(inputFile, rotLib, atLib, pairLib);
    resultList[0]->updateAccuracy(natGraph);

    resultList[0]->print(out);
    cout << "total key num: " << resultList[0]->keyToEnergy.size() << endl;
    clock_t end1 = clock();
	cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;

    out.close();

    delete natGraph;
    delete pairLib;
    delete rotLib;
    delete atLib;
    delete moveLib;
    delete et;
    delete eiLib;
    delete para;
}


