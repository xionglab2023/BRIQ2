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

int testSingleBasePrediction(NuPairMoveSetLibrary* moveLib, EdgeMoveClustersLib* eiLib, RnaEnergyTable* et, BasePairLib* pairLib,RotamerLib* rotLib,AtomLib* atLib, const string& inputFile, const string& outpdb, int pos){
	
    int baseNum = 0;
    double totEne = 0.0;
    double totRms = 0.0;
    NSPtools::InputParser input(inputFile);
    string cst = input.getValue("cst");

    {
        cout << "pos: " << pos << " " << cst[pos] << endl;
        cout << "init graph" << endl;
        NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, eiLib, et);

        cout << "init for single residue prediction" << endl;
        graph->initForSingleResiduePrediction(inputFile, pos);
	    graph->generateRandomEdgePartition(2);
	    SamplingGraph* tree = new SamplingGraph(graph); 
        tree->updatePartitionInfo();
        tree->updateSamplingInfo();

	    //tree->printNodeInfo();
        //tree->printEdgeInfo();
        cout << "run mc" << endl;
	    graphInfo* gi = tree->runAtomicMC();
        baseNum++;
        totEne += gi->ene;
        totRms += gi->rms;
        //graph->printEnergy();
       
        gi->printPDB(outpdb);

        delete gi;
	    delete tree;
	    delete graph;
    }

    return 0;
}

int main(int argc, char** argv){
        

	NuPairMoveSetLibrary* moveLib = new NuPairMoveSetLibrary("xtb", true, 1);
	moveLib->load();
   

    CmdArgs cmdArgs{argc, argv};

    string inputFile = cmdArgs.getValue("-in");
    string outpdb = cmdArgs.getValue("-out");
    int pos = atoi(cmdArgs.getValue("-pos").c_str());

   
    int startID = 0;
	clock_t start = clock();

    ForceFieldPara* para = new ForceFieldPara();
    para->bwTag = "adj";

    para->kStepNum1 = 100;
    para->kStepNum2 = 300;
    para->kStepNum3 = 600;
    para->kNodeFreq = 0.8;

    srand(time(0));


    cout << "load energy table" << endl;
    RnaEnergyTable* et = new RnaEnergyTable(para);
	et->loadAtomicEnergy();

    cout << "init pairLib rotLib atLib" << endl;
    BasePairLib* pairLib = new BasePairLib();
	RotamerLib* rotLib = new RotamerLib();
	AtomLib* atLib = new AtomLib();
    EdgeMoveClustersLib* eiLib = new EdgeMoveClustersLib();

    testSingleBasePrediction(moveLib, eiLib, et, pairLib, rotLib, atLib, inputFile, outpdb, pos);

    delete moveLib;
    delete para;
    delete et;
    delete pairLib;
    delete rotLib;
    delete atLib;
    delete eiLib;
}