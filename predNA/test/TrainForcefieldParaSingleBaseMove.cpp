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
	graph->initRandWeight();
	NuTree* tree = new NuTree(graph);
	graph->MST_kruskal(tree);
	tree->printEdges();
	tree->updateNodeInfo();
	tree->updateEdgeInfo();
	tree->updateSamplingInfo();
	tree->printNodeInfo();    
    return (int)(tree->totalSamp);
}

int testSingleBasePrediction(NuPairMoveSetLibrary* moveLib, RnaEnergyTable* et, BasePairLib* pairLib,RotamerLib* rotLib,AtomLib* atLib, const string& inputFile, double* outEne, double* outRMS, int* outBaseNum, int mpID){
	
    int baseNum = 0;
    double totEne = 0.0;
    double totRms = 0.0;
    NSPtools::InputParser input(inputFile);
    string cst = input.getValue("cst");

    for(int pos=0;pos<cst.length();pos++){
        cout << "pos: " << pos << " " << cst[pos] << endl;
        if(cst[pos] == 'F') continue;
        cout << "init graph" << endl;
        NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, et);

        cout << "init for single residue prediction" << endl;
        graph->initForSingleResiduePrediction(inputFile, pos);
        cout << "init rand weight" << endl;
	    graph->initRandWeight();
        
        cout << "new tree" << endl;
        NuTree* tree = new NuTree(graph);
        cout << "mst" << endl;
	    graph->MST_kruskal(tree);

	    tree->printEdges();
	    tree->updateNodeInfo();
	    tree->updateEdgeInfo();
	    tree->updateSamplingInfo();
	    //tree->printNodeInfo();
        //tree->printEdgeInfo();
        cout << "run mc" << endl;
	    graphInfo* gi = tree->runAtomicMC();
        baseNum++;
        totEne += gi->ene;
        totRms += gi->rms;
       
        delete gi;
	    delete tree;
	    delete graph;
    }

    outEne[mpID] = totEne/baseNum;
    outRMS[mpID] = totRms/baseNum;
    outBaseNum[mpID] = baseNum;

    return 0;
}

int testSingleBasePredictionPrintPDB(NuPairMoveSetLibrary* moveLib, RnaEnergyTable* et, BasePairLib* pairLib,RotamerLib* rotLib,AtomLib* atLib, const string& inputFile, ofstream& out, const string& pdbFile){
	
    NSPtools::InputParser input(inputFile);
    string cst = input.getValue("cst");

    string pdbTag = pdbFile.substr(0, pdbFile.length()-4);
    char xx[200];
    for(int pos=0;pos<cst.length();pos++){
        if(cst[pos] == 'F') continue;
        sprintf(xx, "%s-%d.pdb", pdbTag.c_str(), pos);
        string pdbout = string(xx);
        NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, et);
        graph->initForSingleResiduePrediction(inputFile, pos);
	    graph->initRandWeight();
        NuTree* tree = new NuTree(graph);
	    graph->MST_kruskal(tree);
	    tree->printEdges();
	    tree->updateNodeInfo();
	    tree->updateEdgeInfo();
	    tree->updateSamplingInfo();
        tree->printNodeInfo();
	    graphInfo* gi = tree->runAtomicMC();
        out << pos << " " <<  gi->ene  << " " << gi->rms << endl;
        gi->printPDB(pdbout);
        graph->printEnergy();

        delete gi;
	    delete tree;
	    delete graph;
    }
    return 0;
}

void printHelp(){
    cout << "testSamp -in $INPUTFILE -out $OUTPUTFILE -tag " << endl;
}

int main(int argc, char** argv){
        
    if(argc == 1) {
        printHelp();
        return EXIT_SUCCESS;
    }
    CmdArgs cmdArgs{argc, argv};

    cout << "load move lib: " << endl;
	NuPairMoveSetLibrary* moveLib = new NuPairMoveSetLibrary(true, 1);
	moveLib->load();

    string inputFile = cmdArgs.getValue("-in");
    string tag = cmdArgs.getValue("-tag");
    string outputFile = cmdArgs.getValue("-out");
    string bw = cmdArgs.getValue("-bw");
    //string outpdb = cmdArgs.getValue("-outpdb");

    ofstream out;
    out.open(outputFile.c_str(), ios::out);

    char xx[200];
    
    int mp = 16;
    int startID = 0;
	clock_t start = clock();

    double* outEneList = new double[mp];
    double* outRMSList = new double[mp];
    int* outBaseNum = new int[mp];

    ForceFieldPara* para = new ForceFieldPara();
    para->bwTag = bw;

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

    cout << "task: " << tag << endl;

    cout << "a" << endl;

    if(tag == "clashNb") {
        for(double x= 1.0; x < 10.0;x = x*1.3) {
           
           para->clashLamdaNb = x;

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outEneList[i] = 0.0;
                outRMSList[i] = 0.0;
                outBaseNum[i] = 0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePrediction, moveLib, et, pairLib, rotLib, atLib, inputFile, outEneList, outRMSList, outBaseNum,i-startID);
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
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRms, outBaseNum[0]);
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
                outBaseNum[i] = 0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePrediction, moveLib, et, pairLib, rotLib, atLib, inputFile, outEneList, outRMSList, outBaseNum,i-startID);
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
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRms, outBaseNum[0]);
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
                outBaseNum[i] = 0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePrediction, moveLib, et, pairLib, rotLib, atLib, inputFile, outEneList, outRMSList, outBaseNum,i-startID);
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
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRms, outBaseNum[0]);
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
                outBaseNum[i] = 0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePrediction, moveLib, et, pairLib, rotLib, atLib, inputFile, outEneList, outRMSList, outBaseNum,i-startID);
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
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRms, outBaseNum[0]);
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
                outBaseNum[i] = 0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePrediction, moveLib, et, pairLib, rotLib, atLib, inputFile, outEneList, outRMSList, outBaseNum,i-startID);
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
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRms, outBaseNum[0]);
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
                outBaseNum[i] = 0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePrediction, moveLib, et, pairLib, rotLib, atLib, inputFile, outEneList, outRMSList, outBaseNum,i-startID);
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
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRms, outBaseNum[0]);
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
                outBaseNum[i] = 0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePrediction, moveLib, et, pairLib, rotLib, atLib, inputFile, outEneList, outRMSList, outBaseNum,i-startID);
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
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRms, outBaseNum[0]);
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
                outBaseNum[i] = 0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePrediction, moveLib, et, pairLib, rotLib, atLib, inputFile, outEneList, outRMSList, outBaseNum,i-startID);
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
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRms, outBaseNum[0]);
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
                outBaseNum[i] = 0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePrediction, moveLib, et, pairLib, rotLib, atLib, inputFile, outEneList, outRMSList, outBaseNum,i-startID);
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
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRms, outBaseNum[0]);
            out << string(xx) << endl;

        }

        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "wtBB") {
        
        for(double x= 0.1; x < -2.0;x = x+0.1) {
           
            para->wtRibose = x;
            para->wtPho = x*0.7;

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 
            for(int i=0;i<mp;i++){
                outEneList[i] = 0.0;
                outRMSList[i] = 0.0;
                outBaseNum[i] = 0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePrediction, moveLib, et, pairLib, rotLib, atLib, inputFile, outEneList, outRMSList, outBaseNum,i-startID);
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
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRms, outBaseNum[0]);
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
                outBaseNum[i] = 0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePrediction, moveLib, et, pairLib, rotLib, atLib, inputFile, outEneList, outRMSList, outBaseNum,i-startID);
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
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRms, outBaseNum[0]);
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
                outBaseNum[i] = 0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePrediction, moveLib, et, pairLib, rotLib, atLib, inputFile, outEneList, outRMSList, outBaseNum,i-startID);
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
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRms, outBaseNum[0]);
            out << string(xx) << endl;

        }

        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;       
    }
    else if(tag == "step1"){
    
        for(double x= 50; x < 1000;x = x*1.5) {
           
            para->kStepNum1 = (int)x;
            para->kStepNum2 = (int)((1000-x)*0.5);
            para->kStepNum3 = (int)((1000-x)*0.5);

            cout << "stepNum: " <<  para->kStepNum1 << " " << para->kStepNum2 << " " << para->kStepNum3 << endl;

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outEneList[i] = 0.0;
                outRMSList[i] = 0.0;
                outBaseNum[i] = 0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePrediction, moveLib, et, pairLib, rotLib, atLib, inputFile, outEneList, outRMSList, outBaseNum,i-startID);
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
            
            sprintf(xx, "%8.2f %8.4f %8.4f %d", x, minEne, minRms, outBaseNum[0]);
            out << string(xx) << endl;

        }
    }
    else if(tag == "step2"){
    
        for(double x= 50; x < 1000;x = x*1.5) {
           
            para->kStepNum2 = (int)x;
            para->kStepNum1 = (int)((1000-x)*0.5);
            para->kStepNum3 = (int)((1000-x)*0.5);

            cout << "stepNum: " <<  para->kStepNum1 << " " << para->kStepNum2 << " " << para->kStepNum3 << endl;


            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outEneList[i] = 0.0;
                outRMSList[i] = 0.0;
                outBaseNum[i] = 0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePrediction, moveLib, et, pairLib, rotLib, atLib, inputFile, outEneList, outRMSList, outBaseNum,i-startID);
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
            
           sprintf(xx, "%8.2f %8.4f %8.4f %d", x, minEne, minRms, outBaseNum[0]);
            out << string(xx) << endl;

        }

        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;       
    }
    else if(tag == "step3"){
    
        for(double x= 50; x < 1000;x = x*1.5) {
           
            para->kStepNum3 = (int)x;
            para->kStepNum1 = (int)((1000-x)*0.5);
            para->kStepNum2 = (int)((1000-x)*0.5);

            cout << "stepNum: " <<  para->kStepNum1 << " " << para->kStepNum2 << " " << para->kStepNum3 << endl;

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outEneList[i] = 0.0;
                outRMSList[i] = 0.0;
                outBaseNum[i] = 0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePrediction, moveLib, et, pairLib, rotLib, atLib, inputFile, outEneList, outRMSList, outBaseNum,i-startID);
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
            
            sprintf(xx, "%8.2f %8.4f %8.4f %d", x, minEne, minRms, outBaseNum[0]);
            out << string(xx) << endl;

        }
   }
    else if(tag == "kNode"){
    
        for(double x= 0.1; x < 5.0;x = x*1.5) {
           
            para->kStepNum3 = 333;
            para->kStepNum1 = 333;
            para->kStepNum2 = 333;

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outEneList[i] = 0.0;
                outRMSList[i] = 0.0;
                outBaseNum[i] = 0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePrediction, moveLib, et, pairLib, rotLib, atLib, inputFile, outEneList, outRMSList, outBaseNum,i-startID);
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
            
            sprintf(xx, "%8.2f %8.4f %8.4f %d", x, minEne, minRms, outBaseNum[0]);
            out << string(xx) << endl;

        }
    }
    else if(tag == "test") {
       testSingleBasePrediction(moveLib, et, pairLib, rotLib, atLib, inputFile, outEneList, outRMSList, outBaseNum, 0);
    }
    else if(tag == "bw") {

        {
         
            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outEneList[i] = 0.0;
                outRMSList[i] = 0.0;
                outBaseNum[i] = 0;
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePrediction, moveLib, et, pairLib, rotLib, atLib, inputFile, outEneList, outRMSList, outBaseNum,i-startID);
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
            
            sprintf(xx, "bw %8.4f %8.4f %d",  minEne, minRms, outBaseNum[0]);
            out << string(xx) << endl;
        }
    }
    else if(tag == "pdb") {
        //testSingleBasePredictionPrintPDB(moveLib, et, pairLib, rotLib, atLib, inputFile, out, outpdb);
    }

    out.close();
    delete pairLib;
    delete rotLib;
    delete atLib;
    delete moveLib;
    delete et;
    delete para;

}
