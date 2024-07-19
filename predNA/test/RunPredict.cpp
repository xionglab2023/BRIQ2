/*
 * RunPredict.cpp
 *
 */

#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include "model/StructureModel.h"
#include "predNA/NuGraph.h"
#include "predNA/NuSampling.h"
#include "tools/CmdArgs.h"
#include "tools/ThreadPool.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredNA;
using namespace std;
using namespace NSPtools;
using namespace NSPthread;

int runRefinement(NuPairMoveSetLibrary* moveLib, EdgeInformationLib* eiLib, RnaEnergyTable* et, const string& inputFile, const string& outFile, int randSeed, double kStep){

    srand(randSeed);

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

	BasePairLib* pairLib = new BasePairLib("stat");
	RotamerLib* rotLib = new RotamerLib();
	AtomLib* atLib = new AtomLib();
	NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, eiLib, et);
    graph->initForMC(inputFile);
	graph->initRandWeight();
	NuTree* tree = new NuTree(graph);
	graph->MST_kruskal(tree);
	tree->updateNodeInfo(1.0, 1.0);
	tree->updateEdgeInfo(1.0, 1.0);
	tree->updateSamplingInfo();

	cout << "run MC" << endl;
	graphInfo* gi = tree->runAtomicMC();
    gi->printPDB(outFile);

    delete gi;
    delete pairLib;
    delete rotLib;
    delete atLib;
	delete tree;
	delete graph;
    return 0;
}

int runPredict(NuPairMoveSetLibrary* moveLib, EdgeInformationLib* eiLib, RnaEnergyTable* et, const string& inputFile, const string& outFile, int randSeed, double kStep){
    srand(randSeed);

	et->para->kStepNum1 = (int)(et->para->kStepNum1*kStep);
	et->para->kStepNum2 = (int)(et->para->kStepNum2*kStep);
	et->para->kStepNum3 = (int)(et->para->kStepNum3*kStep);
	
	BasePairLib* pairLib = new BasePairLib("stat");
	RotamerLib* rotLib = new RotamerLib();
	AtomLib* atLib = new AtomLib();
	NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, eiLib, et);
    graph->initForMC(inputFile);
	graph->initRandWeight();
	NuTree* tree = new NuTree(graph);
	graph->MST_kruskal(tree);
	tree->updateNodeInfo(1.0, 1.0);
	tree->updateEdgeInfo(1.0, 1.0);
	tree->updateSamplingInfo();

	cout << "run MC" << endl;
	graphInfo* gi = tree->runAtomicMC();
    gi->printPDB(outFile);

    delete gi;
    delete pairLib;
    delete rotLib;
    delete atLib;
	delete tree;
	delete graph;

    return 0;
}

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

int runCGMC(NuPairMoveSetLibrary* moveLib, RnaEnergyTable* et,EdgeInformationLib* eiLib, const string& inputFile, cgResult* result, int modelNum, int randSeed, double kStep){

 	srand(randSeed);

	et->para->kStepNum1CG = (int)(et->para->kStepNum1CG*kStep);
	et->para->kStepNum2CG = (int)(et->para->kStepNum2CG*kStep);
	et->para->kStepNum3CG = (int)(et->para->kStepNum3CG*kStep);

    BasePairLib* pairLib = new BasePairLib("stat");
	RotamerLib* rotLib = new RotamerLib(et->para);
    AtomLib* atLib = new AtomLib();


    cout << "init graph" << endl;
	NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, eiLib, et);

	graph->initForCGMC(inputFile);
	graph->initRandWeight();
	NuTree* tree = new NuTree(graph);
	graph->MST_kruskal(tree);
	tree->updateNodeInfoCG(1.0, 1.0);
	tree->updateEdgeInfoCG(1.0, 1.0);
	tree->updateSamplingInfo();

	clock_t start = clock();
    cout << "run mc" << endl;

    NuSampling* samp = new NuSampling(graph, tree);

    samp->runCoarseGrainedMC(result->keyToEnergy, modelNum);

	delete pairLib;
    delete rotLib;
    delete atLib;
	delete tree;
	delete graph;
    return 0;
}

int cgKeyToAllAtomPDB(NuPairMoveSetLibrary* moveLib, EdgeInformationLib* eiLib, RnaEnergyTable* et, const string& inputFile, const string& outFile, int randSeed, double kStep){
    srand(randSeed);

   	et->para->T0 = 10.0;
    et->para->T1 = 1.0;
    et->para->T2 = 0.1;
   	et->para->T3 = 0.01;
    et->para->clashRescale = 0.2;
    et->para->connectRescale = 0.5;
	et->para->kStepNum1 = 200;
	et->para->kStepNum2 = 100;
	et->para->kStepNum3 = 100;
	et->para->withRandomInit = false;

	et->para->kStepNum1 = (int)(et->para->kStepNum1*kStep);
	et->para->kStepNum2 = (int)(et->para->kStepNum2*kStep);
	et->para->kStepNum3 = (int)(et->para->kStepNum3*kStep);
	

	BasePairLib* pairLib = new BasePairLib("stat");
	RotamerLib* rotLib = new RotamerLib();
	AtomLib* atLib = new AtomLib();
	NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, eiLib, et);
	cout << "init for mc" << endl;
    graph->initForMC(inputFile);
	cout << "init rand weight" << endl;
	graph->initRandWeight();
	graph->printAllEdge();
	cout << "init tree" << endl;
	NuTree* tree = new NuTree(graph);
	cout << "MST" << endl;
	graph->MST_kruskal(tree);
	cout << "update node info" << endl;
	tree->updateNodeInfo(1.0, 1.0);
	cout << "update edge info" << endl;
	tree->updateEdgeInfo(1.0, 1.0);
	cout << "update sampling info" << endl;
	tree->updateSamplingInfo();

	cout << "run MC" << endl;
	tree->runAtomicMC();
	
	graph->initNearestNativeEdge();
	graph->initRandWeight();
	graph->MST_kruskal(tree);
	tree->updateNodeInfo(1.0, 1.0);
	tree->updateEdgeInfo(1.0, 1.0);
	tree->updateSamplingInfo();

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

	et->para->kStepNum1 = (int)(et->para->kStepNum1*kStep);
	et->para->kStepNum2 = (int)(et->para->kStepNum2*kStep);
	et->para->kStepNum3 = (int)(et->para->kStepNum3*kStep);

	cout << "run MC2" << endl;
	graphInfo* gi = tree->runAtomicMC();

    gi->printPDB(outFile);

    delete gi;
    delete pairLib;
    delete rotLib;
    delete atLib;
	delete tree;
	delete graph;
    return 0;
}

void printHelp(){
	cout << "Usage: " << endl;
	cout << "run refinement:" <<endl;
	cout << "briqx_run -in $INPUTFILE -out $OUTPUT -seed $RANDOM_SEED" << endl;
	cout << "run predict: " << endl;
	cout << "briqx_run -in $INPUTFILE -out $OUTPUT -seed $RANDOM_SEED" << endl;
	cout << "run cgModeling: " << endl;
	cout << "briqx_run -in $INPUTFILE -out $OUTPUT -seed $RANDOM_SEED -n $KEYNum -mp $THREADNum" << endl;
}

int main(int argc, char** argv){

    if(argc == 1) {
        printHelp();
        return EXIT_SUCCESS;
    }

	clock_t start = clock();

    CmdArgs cmdArgs{argc, argv};
	if(cmdArgs.specifiedOption("-h")){
		printHelp();
        return EXIT_SUCCESS;
	}

	string inputFile = cmdArgs.getValue("-in");
    string outputFile = cmdArgs.getValue("-out");
	int seed = 0;
	if(cmdArgs.specifiedOption("-seed")) {
		seed = atoi(cmdArgs.getValue("-seed").c_str());
	}
	else {
		seed = time(0);
	}

	double kStep = 1.0;
	if(cmdArgs.specifiedOption("-kStep")){
		kStep = atof(cmdArgs.getValue("-kStep").c_str());
	}

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
	NSPtools::InputParser input(inputFile);
	
	string task = input.getValue("task");
	string key = input.getValue("key");

	if(task == "refinement"){
		et->loadAtomicEnergy();
		runRefinement(moveLib, eiLib, et, inputFile, outputFile, seed, kStep);
	}
	else if(task == "predict" && key.length() == 0) {
		et->loadAtomicEnergy();
		runPredict(moveLib, eiLib, et, inputFile, outputFile, seed, kStep);
	}
	else if(task == "predict"){
		et->loadAtomicEnergy();
		cgKeyToAllAtomPDB(moveLib, eiLib, et, inputFile, outputFile, seed, kStep);
	}
	else if(task == "cgModeling"){
		
		et->loadCoarseGrainedEnergy();
		int modelNum = 10;
		int mp = 1;
		if(cmdArgs.specifiedOption("-n"))
			modelNum = atoi(cmdArgs.getValue("-n").c_str());
		else if(cmdArgs.specifiedOption("-mp"))
			mp = atoi(cmdArgs.getValue("-mp").c_str());
    	

    	ofstream out;
   	 	out.open(outputFile.c_str(), ios::out);

    	int startID = 0;
    

		

    	shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
    	size_t jid = 0;
    	char xx[200];

    	cgResult* resultList[mp];
    	for(int i=0;i<mp;i++){
        	resultList[i] = new cgResult();
    	}

    	for(int i=startID;i<startID+mp;i++) {
        	shared_ptr<IntFuncTask> request(new IntFuncTask);

        	request->asynBind(runCGMC, moveLib, et, eiLib, inputFile,  resultList[i-startID],modelNum, time(0)+i, kStep);
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

    	resultList[0]->print(out);
    	out.close();
	}

	delete pairLib;
	delete rotLib;
	delete atLib;
	delete moveLib;
	delete eiLib;
	delete para;
	delete et;

	clock_t end1 = clock();
	cout << "time1: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
}


