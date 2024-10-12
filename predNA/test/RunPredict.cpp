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

int runRefinement(NuPairMoveSetLibrary* moveLib, EdgeInformationLib* eiLib, RnaEnergyTable* et, const string& inputFile, const string& outFile, int randSeed, double kStep, bool printDetail){

    srand(randSeed);

   	et->para->T0 = 1.0;
    et->para->T1 = 0.2;
    et->para->T2 = 0.05;
   	et->para->T3 = 0.01;
    et->para->clashRescale = 0.2;
    et->para->connectRescale = 0.5;
	et->para->kStepNum1 = (int)(200*kStep);
	et->para->kStepNum2 = (int)(100*kStep);
	et->para->kStepNum3 = (int)(100*kStep);
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
	if(printDetail){
		cout << "detail energy: " << endl;
		graph->printEnergy();
	}

    delete gi;
    delete pairLib;
    delete rotLib;
    delete atLib;
	delete tree;
	delete graph;
    return 0;
}

int runRefinementFromPDB(NuPairMoveSetLibrary* moveLib, EdgeInformationLib* eiLib, RnaEnergyTable* et, const string& pdbFile, const string& outFile, int randSeed, double kStep, const string& cnt){

    srand(randSeed);

   	et->para->T0 = 0.3;
    et->para->T1 = 0.07;
    et->para->T2 = 0.015;
   	et->para->T3 = 0.002;
    et->para->clashRescale = 0.2;
    et->para->connectRescale = 0.4;
	et->para->kStepNum1 = (int)(500*kStep);
	et->para->kStepNum2 = (int)(500*kStep);
	et->para->kStepNum3 = (int)(500*kStep);
	et->para->withRandomInit = false;

	cout << "init libs" << endl;

	BasePairLib* pairLib = new BasePairLib("stat");
	RotamerLib* rotLib = new RotamerLib();
	AtomLib* atLib = new AtomLib();

	cout << "init pdb" << endl;
	cout << "pdbFile: " << pdbFile << endl;
	RNAPDB* pdb = new RNAPDB(pdbFile);
	cout << "init graph" <<endl;
	NuGraph* graph = new NuGraph(pdb, rotLib, atLib, pairLib, moveLib, eiLib, et, cnt);
	cout << "init weight" << endl;
	graph->initRandWeight();
	NuTree* tree = new NuTree(graph);
	graph->MST_kruskal(tree);
	tree->updateNodeInfo(1.0, 1.0);
	tree->updateEdgeInfo(1.0, 1.0);
	tree->updateSamplingInfo();

	tree->printNodeInfo();
	tree->printEdgeInfo();

	cout << "run MC" << endl;
	graphInfo* gi = tree->runAtomicMC();
    gi->printPDB(outFile);
	graph->printBaseEnergyList(outFile);
	graph->printPairwiseEnergy(outFile);

	delete pdb;
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

int runCGMC(NuPairMoveSetLibrary* moveLib, RnaEnergyTable* et,EdgeInformationLib* eiLib, const string& inputFile, const string& outFile, int modelNum, int randSeed, double kStep){

 	srand(randSeed);

	et->para->kStepNum1CG = (int)(et->para->kStepNum1CG*kStep);
	et->para->kStepNum2CG = (int)(et->para->kStepNum2CG*kStep);
	et->para->kStepNum3CG = (int)(et->para->kStepNum3CG*kStep);

    BasePairLib* pairLib = new BasePairLib("stat");
	RotamerLib* rotLib = new RotamerLib(et->para);
    AtomLib* atLib = new AtomLib();

	NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, eiLib, et);

	graph->initForCGMC(inputFile);
	graph->initRandWeight();
	NuTree* tree = new NuTree(graph);
	graph->MST_kruskal(tree);
	tree->updateNodeInfoCG(1.0, 1.0);
	tree->updateEdgeInfoCG(1.0, 1.0);
	tree->updateSamplingInfo();

	clock_t start = clock();

	
    NuSampling* samp = new NuSampling(graph, tree);
	cout << "runMC" << endl;
    samp->runCoarseGrainedMC(outFile, modelNum);

	delete pairLib;
    delete rotLib;
    delete atLib;
	delete tree;
	delete graph;
    return 0;
}

int cgKeyToAllAtomPDB(NuPairMoveSetLibrary* moveLib, EdgeInformationLib* eiLib, RnaEnergyTable* et, const string& inputFile, const string& outFile, int modelNum, int randSeed, double kStep){
    srand(randSeed);
	BasePairLib* pairLib = new BasePairLib("stat");
	RotamerLib* rotLib = new RotamerLib();
	AtomLib* atLib = new AtomLib();
	NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, eiLib, et);
	graph->initForMC(inputFile);
	NuTree* tree = new NuTree(graph);

	ofstream out;
	out.open(outFile.c_str(), ios::out);
	if(!out.is_open()) {
		cout << "fail to open file: " << outFile << endl;
		exit(0);
	}
	out.close();

	for(int i=0;i<modelNum;i++){

		et->para->T0 = 10.0;
    	et->para->T1 = 1.0;
    	et->para->T2 = 0.1;
   		et->para->T3 = 0.01;
    	et->para->clashRescale = 0.2;
    	et->para->connectRescale = 0.5;
		et->para->kStepNum1 = 200;
		et->para->kStepNum2 = 100;
		et->para->kStepNum3 = 100;
		et->para->withRandomInit = true;

		et->para->kStepNum1 = (int)(et->para->kStepNum1*kStep);
		et->para->kStepNum2 = (int)(et->para->kStepNum2*kStep);
		et->para->kStepNum3 = (int)(et->para->kStepNum3*kStep);

		graph->initRandWeight();
		graph->printAllEdge();

		graph->MST_kruskal(tree);
		tree->updateNodeInfo(1.0, 1.0);
		tree->updateEdgeInfo(1.0, 1.0);
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

		tree->printTreeInfomation(outFile);
		delete gi;
	}


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
	cout << "briqx_run -pdb $PDBFILE -out OUTPDB -seed $RANDOM_SEED" << endl;
	cout << "briqx_run -pdb $PDBFILE -out OUTPDB -seed $RANDOM_SEED -cnt \"X$CNT\"" << endl;  
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

	string inputFile = "";
	string pdbFile = "";
	if(cmdArgs.specifiedOption("-in")){
		inputFile = cmdArgs.getValue("-in");
	}
	if(cmdArgs.specifiedOption("-pdb")){
		pdbFile = cmdArgs.getValue("-pdb");
	}    
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

	if(cmdArgs.specifiedOption("-pdb")){
		et->loadAtomicEnergy();
		cout << "run refinement from PDB" << endl;
		string cnt = "";
		if(cmdArgs.specifiedOption("-cnt"))
			cnt = cmdArgs.getValue("-cnt");

		if(cnt.length() > 1)
			cnt = cnt.substr(1, cnt.length()-1);
			
		cout << "cnt: " << cnt << endl;
		runRefinementFromPDB(moveLib, eiLib, et, pdbFile, outputFile, seed, kStep, cnt);
	}
	else {
		if(inputFile.length() == 0){
			cout << "please specify pdb file or inputFile" << endl;
			exit(0);
		}
		NSPtools::InputParser input(inputFile);
	
		string task = input.getValue("task");
		string key = input.getValue("key");

		bool printEnergyDetail = false;
		if(cmdArgs.specifiedOption("-v")){
			printEnergyDetail = true;
		}

		if(task == "refinement"){
			et->loadAtomicEnergy();
			runRefinement(moveLib, eiLib, et, inputFile, outputFile, seed, kStep, printEnergyDetail);
		}
		else if(task == "predict" && key.length() == 0) {
			et->loadAtomicEnergy();
			runPredict(moveLib, eiLib, et, inputFile, outputFile, seed, kStep);
		}
		else if(task == "predict"){
			int modelNum = 10;
			if(cmdArgs.specifiedOption("-n"))
				modelNum = atoi(cmdArgs.getValue("-n").c_str());

			et->loadAtomicEnergy();
			cgKeyToAllAtomPDB(moveLib, eiLib, et, inputFile, outputFile, modelNum, seed, kStep);
		}
		else if(task == "cgModeling"){

			int modelNum = 10;
			int mp = 1;
			if(cmdArgs.specifiedOption("-n"))
				modelNum = atoi(cmdArgs.getValue("-n").c_str());

			et->loadCoarseGrainedEnergy();
			runCGMC(moveLib, et, eiLib, inputFile, outputFile, modelNum, seed, kStep);
		}

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


