
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


int testSamplingStepNum(NuPairMoveSetLibrary* moveLib, RnaEnergyTable* et, const string& inputFile){
    BasePairLib* pairLib = new BasePairLib();
	RotamerLib* rotLib = new RotamerLib();
	AtomLib* atLib = new AtomLib();
	NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, et);
    graph->initForMC(inputFile);
	graph->initRandWeight();
	NuTree* tree = new NuTree(graph);
	graph->MST_kruskal(tree);
	tree->printEdges();
	tree->updateNodeInfo();
	tree->updateEdgeInfo();
	tree->updateSamplingInfo();
	tree->printNodeInfo();    
    return (int)tree->totalSamp;
}

int testRefinement(NuPairMoveSetLibrary* moveLib, RnaEnergyTable* et, const string& inputFile, double* outEne, double* outRMS, int mpID){

    BasePairLib* pairLib = new BasePairLib();
	RotamerLib* rotLib = new RotamerLib();
	AtomLib* atLib = new AtomLib();
	NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, et);
    graph->initForMC(inputFile);
	graph->initRandWeight();
	NuTree* tree = new NuTree(graph);
	graph->MST_kruskal(tree);
	tree->printEdges();
	tree->updateNodeInfo();
	tree->updateEdgeInfo();
	tree->updateSamplingInfo();
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

	NuPairMoveSetLibrary* moveLib = new NuPairMoveSetLibrary(true, 1);
	moveLib->load();

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
    para->kStepNum1 = 600;
    para->kStepNum2 = 600;
    para->kStepNum3 = 600;
    RnaEnergyTable* et = new RnaEnergyTable(para);
	et->loadAtomicEnergy();

    if(tag == "clashNb") {
        for(double x= 1.0; x < 10.0;x = x*1.3) {
           
           para->clashLamdaNb = x;

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outEneList[i] = 0.0;
                outRMSList[i] = 0.0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, et, inputFile, outEneList, outRMSList, i-startID);
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

            sprintf(xx, "%4.2f %8.4f %8.4f", para->clashLamdaNb, minEne, minRms);
            out << string(xx) << endl;
        }

        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "clashNnb") {
        for(double x= 1.0; x < 10.0;x = x*1.3) {
           
           para->clashLamdaNnb = x;

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outEneList[i] = 0.0;
                outRMSList[i] = 0.0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, et, inputFile, outEneList, outRMSList, i-startID);
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
            sprintf(xx, "%4.2f %8.4f %8.4f", para->clashLamdaNnb, minEne, minRms);
            out << string(xx) << endl;
        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "hbLamda1") {
        for(double x= 0.03; x < 0.5;x = x*1.3) {
           
           para->hbLamda1 = x;

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outEneList[i] = 0.0;
                outRMSList[i] = 0.0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, et, inputFile, outEneList, outRMSList, i-startID);
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
            sprintf(xx, "%5.3f %8.4f %8.4f", para->hbLamda1, minEne, minRms);
            out << string(xx) << endl;
        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "hbLamda2") {
        for(double x= 0.1; x < 1.0;x = x*1.3) {
           
           para->hbLamda2 = x;

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outEneList[i] = 0.0;
                outRMSList[i] = 0.0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, et, inputFile, outEneList, outRMSList, i-startID);
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
            sprintf(xx, "%5.3f %8.4f %8.4f", x, minEne, minRms);
            out << string(xx) << endl;
        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "wtHb") {
        for(double x= 0.3; x < 5.0;x = x*1.3) {
           
           para->wtHb = x;

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outEneList[i] = 0.0;
                outRMSList[i] = 0.0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, et, inputFile, outEneList, outRMSList, i-startID);
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
            sprintf(xx, "%5.3f %8.4f %8.4f", x, minEne, minRms);
            out << string(xx) << endl;
        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "wtBp1") {
        for(double x= 0.3; x < 6.0;x = x*1.3) {
           
           para->wtBp1 = x;

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outEneList[i] = 0.0;
                outRMSList[i] = 0.0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, et, inputFile, outEneList, outRMSList, i-startID);
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
            sprintf(xx, "%5.3f %8.4f %8.4f", x, minEne, minRms);
            out << string(xx) << endl;
        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "wtBp2") {
        for(double x= 0.3; x < 6.0;x = x*1.3) {
           
           para->wtBp2 = x;

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outEneList[i] = 0.0;
                outRMSList[i] = 0.0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, et, inputFile, outEneList, outRMSList, i-startID);
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
            sprintf(xx, "%5.3f %8.4f %8.4f", x, minEne, minRms);
            out << string(xx) << endl;
        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "wtO4O2C2Nb") {
        for(double x= 0.0; x < 2.0;x = x+0.2) {
           
           para->wtO4O2C2Nb = x;

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outEneList[i] = 0.0;
                outRMSList[i] = 0.0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, et, inputFile, outEneList, outRMSList, i-startID);
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
            sprintf(xx, "%5.3f %8.4f %8.4f", x, minEne, minRms);
            out << string(xx) << endl;
        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "wtO4O2C2Nnb") {
        for(double x= 0.0; x < 1.21;x = x+0.2) {
           
           para->wtO4O2C2Nnb = x;

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outEneList[i] = 0.0;
                outRMSList[i] = 0.0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, et, inputFile, outEneList, outRMSList, i-startID);
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
            sprintf(xx, "%5.3f %8.4f %8.4f", x, minEne, minRms);
            out << string(xx) << endl;
        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "wtBB") {
        for(double x= 0.1; x < 2.0;x = x+0.1) {
           
            para->wtRibose = x;
            para->wtPho = x*0.7;

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outEneList[i] = 0.0;
                outRMSList[i] = 0.0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, et, inputFile, outEneList, outRMSList, i-startID);
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
            sprintf(xx, "%5.3f %8.4f %8.4f", x, minEne, minRms);
            out << string(xx) << endl;
        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "wtOxy") {
        for(double x= 0.1; x < 2.0;x = x+0.1) {
            for(int k=0;k<16;k++) {
                para->wtRiboseOxy[k] = x;
            }

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outEneList[i] = 0.0;
                outRMSList[i] = 0.0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, et, inputFile, outEneList, outRMSList, i-startID);
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
            sprintf(xx, "%5.3f %8.4f %8.4f", x, minEne, minRms);
            out << string(xx) << endl;
        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "T0"){
        for(double x= 2.0; x < 50.0;x = x*1.3) {
            para->T0 = x;
            para->T1 = para->T0*0.1;
            para->T2 = para->T1*0.1;
            para->T3 = para->T2*0.1;

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outEneList[i] = 0.0;
                outRMSList[i] = 0.0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, et, inputFile, outEneList, outRMSList, i-startID);
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
            sprintf(xx, "%5.3f %8.4f %8.4f", x, minEne, minRms);
            out << string(xx) << endl;
        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;        
    }

    else if(tag == "test") {
        int stepNum = testSamplingStepNum(moveLib, et, inputFile);
        out << stepNum << endl;
    }

    out.close();
    delete moveLib;
    delete et;
    delete para;

}