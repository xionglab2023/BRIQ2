
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


int testRefinement(NuPairMoveSetLibrary* moveLib, EdgeMoveClustersLib* eiLib, RnaEnergyTable* et, const string& inputFile, double* outEne, double* outRMS, int mpID){

    BasePairLib* pairLib = new BasePairLib();
	RotamerLib* rotLib = new RotamerLib();
	AtomLib* atLib = new AtomLib();
	NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, eiLib, et);
	graph->generateRandomEdgePartition(2);
	SamplingGraph* tree = new SamplingGraph(graph);
    tree->updatePartitionInfo();
    tree->updateSamplingInfo();

	for(int i=0;i<graph->seqLen;i++){
		graph->allNodes[i]->printNodeInfo();
	}


	for(int i=0;i<tree->geList.size();i++){
		cout << "edge: " << tree->geList[i]->indexA << "-" << tree->geList[i]->indexB << endl;
		tree->geList[i]->printPartition();
	}

	tree->printNodeInfo();

	cout << "run MC" << endl;
	graphInfo* gi = tree->runAtomicMC();
    outRMS[mpID] = gi->rms;
    outEne[mpID] = gi->ene;

    delete gi;

    delete pairLib;
    delete rotLib;
    delete atLib;
	delete tree;
	delete graph;

    return 0;
}

void printHelp(){
    cout << "testSamp -in $INPUTFILE -out $OUTPUTFILE -mp 64" << endl;
}

int main(int argc, char** argv){
        
    if(argc == 1) {
        printHelp();
        return EXIT_SUCCESS;
    }
    CmdArgs cmdArgs{argc, argv};

	NuPairMoveSetLibrary* moveLib = new NuPairMoveSetLibrary("xtb", true, 1);
	moveLib->load();

    BasePairLib* pairLib = new BasePairLib();
    EdgeMoveClustersLib* eiLib = new EdgeMoveClustersLib();

    string inputFile = cmdArgs.getValue("-in");
    string tag = cmdArgs.getValue("-tag");
    string outputFile = cmdArgs.getValue("-out");

    ofstream out;
    out.open(outputFile.c_str(), ios::out);

    char xx[200];
    
    int mp = 16;
    int startID = 0;
	clock_t start = clock();

    double* outEneList = new double[mp];
    double* outRMSList = new double[mp];

    ForceFieldPara* para = new ForceFieldPara();
    RnaEnergyTable* et = new RnaEnergyTable(para);
	et->loadAtomicEnergy();

    if(tag == "stepNum") {
        for(double n= 10.0; n < 10250;n = n*2.0) {
            para->T0 = 15.0;
            para->T1 = 1.5;
            para->T2 = 0.15;
            para->T3 = 0.015;
            para->kStepNum1 = (int)n;
            para->kStepNum2 = (int)n;
            para->kStepNum3 = (int)n;

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outEneList[i] = 0.0;
                outRMSList[i] = 0.0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, eiLib, et, inputFile, outEneList, outRMSList, i-startID);
                jid++;
                thrPool->addTask(request);
            }

            while(true) {
                sleep(1);
                if(thrPool->getTaskCount() == 0) {
                    break;
                }
             }

            double minEne = 9999.9;
            for(int i=0;i<mp;i++){
                if(outEneList[i] < minEne){
                    minEne = outEneList[i];
                }
            }

            double minRms = 9999.9;
            for(int i=0;i<mp;i++){
                if(outRMSList[i] < minRms){
                    minRms = outRMSList[i];
                }   
            }

            sprintf(xx, "%d %8.4f %8.4f", para->kStepNum1, minEne, minRms);
            out << string(xx) << endl;
        }

        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "T2") {

        //for(int k=2000;k<100000;k=k*2)
        for(double T0=5.0;T0<30.0;T0=T0*1.3){
            double T1 = T0*0.1;
            double T2 = T1*0.1;
            double T3 = T2*0.1;
            int n = 0;
            for(double T=T0;T>T1;T=T*0.95){
                n++;
            }
            for(double T=T1;T>T2;T=T*0.95){
                n++;
            }
            for(double T=T2;T>T3;T=T*0.95){
                n++;
            }

            para->kStepNum1 = 20000/n;
            para->kStepNum2 = 20000/n;
            para->kStepNum3 = 20000/n;

            para->T0 = T0;
            para->T1 = T1;
            para->T2 = T2;
            para->T3 = T3;

   

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outEneList[i] = 0.0;
                outRMSList[i] = 0.0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, eiLib, et, inputFile, outEneList, outRMSList, i-startID);
                jid++;
                thrPool->addTask(request);
            }

            while(true) {
                sleep(1);
                if(thrPool->getTaskCount() == 0) {
                    break;
                }
            }

            double minEne = 9999.9;
            for(int i=0;i<mp;i++){
                if(outEneList[i] < minEne){
                    minEne = outEneList[i];
                }   
            }

            double minRms = 9999.9;
            for(int i=0;i<mp;i++){
                if(outRMSList[i] < minRms){
                    minRms = outRMSList[i];
                }   
            }

            sprintf(xx, "%8.6f %8.6f %8.6f %8.6f %6d %8.4f %8.4f", T0, T1, T2, T3, para->kStepNum1, minEne, minRms);
            out << string(xx) << endl;             

        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }

    out.close();
    delete moveLib;
    delete et;
    delete para;

}
