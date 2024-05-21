
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

class motifPredictionResult{
public:
    double rms;
    double ene;

    motifPredictionResult(){
        
    }

    void setValue( double ene, double rms){
       this->ene = ene;
       this->rms = rms;
    }

    void mergeResult(motifPredictionResult* other){
        if(other->ene < this->ene){
            this->rms = other->rms;
            this->ene = other->ene;
        }
    }

};


int testSamplingStepNum(NuPairMoveSetLibrary* moveLib, EdgeInformationLib* eiLib, RnaEnergyTable* et, const string& inputFile){
    BasePairLib* pairLib = new BasePairLib();
	RotamerLib* rotLib = new RotamerLib();
	AtomLib* atLib = new AtomLib();
	NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, eiLib, et);
    graph->initForMC(inputFile);
	graph->initRandWeight();
	NuTree* tree = new NuTree(graph);
	graph->MST_kruskal(tree);
	tree->printEdges();
	tree->updateNodeInfo(1.0, 1.0);
	tree->updateEdgeInfo(1.0, 1.0);
	tree->updateSamplingInfo();
	tree->printNodeInfo();    
    return (int)tree->totalSamp;
}

int testRefinement(NuPairMoveSetLibrary* moveLib, EdgeInformationLib* eiLib, RnaEnergyTable* et, const string& inputFile, motifPredictionResult** results, int mpID){

    BasePairLib* pairLib = new BasePairLib("stat");
	RotamerLib* rotLib = new RotamerLib();
	AtomLib* atLib = new AtomLib();
	NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, eiLib, et);
    graph->initForMC(inputFile);
	graph->initRandWeight();
	NuTree* tree = new NuTree(graph);
	graph->MST_kruskal(tree);
	tree->printEdges();
	tree->updateNodeInfo(1.0, 1.0);
	tree->updateEdgeInfo(1.0, 1.0);
	tree->updateSamplingInfo();
	tree->printNodeInfo();

	graphInfo* gi = tree->runAtomicMC();

    results[mpID]->setValue(gi->ene, gi->rms);

    delete gi;
    delete pairLib;
    delete rotLib;
    delete atLib;
	delete tree;
	delete graph;

    return 0;
}

void printHelp(){
    cout << "testSamp -in $INPUTFILE -out $OUTPUTFILE -mp 16" << endl;
}

int main(int argc, char** argv){
        
    if(argc == 1) {
        printHelp();
        return EXIT_SUCCESS;
    }
    CmdArgs cmdArgs{argc, argv};

	NuPairMoveSetLibrary* moveLib = new NuPairMoveSetLibrary("stat", true, 1);
	moveLib->load();
    EdgeInformationLib* eiLib = new EdgeInformationLib();

    string inputFile = cmdArgs.getValue("-in");
    string tag = cmdArgs.getValue("-tag");
    string outputFile = cmdArgs.getValue("-out");

    double y = atof(cmdArgs.getValue("-x").c_str());

    ofstream out;
    out.open(outputFile.c_str(), ios::out);

    char xx[200];
    
    int mp = 1;
    int startID = 0;
	clock_t start = clock();

    motifPredictionResult** results = new motifPredictionResult*[mp];
    for(int i=0;i<mp;i++){
        results[i] = new motifPredictionResult();
    }

    ForceFieldPara* para = new ForceFieldPara();
    para->kStepNum1 = 600;
    para->kStepNum2 = 600;
    para->kStepNum3 = 600;

    para->libType = "stat";

    RnaEnergyTable* et = new RnaEnergyTable(para);
	et->loadAtomicEnergy();

    if(tag == "clashNb") {
        for(double x= 1.0; x < 10.0;x = x*1.3) {
           
           para->clashLamdaNb = x;
           et->updateAtomic(para);

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 
  
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, eiLib, et, inputFile, results, i-startID);
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
                results[0]->mergeResult(results[i]);
            }

            double minEne = results[0]->ene;
            double minRms = results[0]->rms;


            sprintf(xx, "%4.2f %8.4f %8.4f", para->clashLamdaNb, minEne, minRms);
            out << string(xx) << endl;
        }

        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "clashNnb") {
        for(double x= 1.0; x < 10.0;x = x*1.3) {
           
           para->clashLamdaNnb = x;
            et->updateAtomic(para);
            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 
 
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, eiLib, et, inputFile, results, i-startID);
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
                results[0]->mergeResult(results[i]);
            }

            double minEne = results[0]->ene;
            double minRms = results[0]->rms;
            sprintf(xx, "%4.2f %8.4f %8.4f", para->clashLamdaNnb, minEne, minRms);
            out << string(xx) << endl;
        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "hbLamda1") {
        for(double x= 0.03; x < 0.5;x = x*1.3) {
           
           para->hbLamda1 = x;
           et->updateAtomic(para);

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, eiLib, et, inputFile, results, i-startID);
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
                results[0]->mergeResult(results[i]);
            }

            double minEne = results[0]->ene;
            double minRms = results[0]->rms;
            sprintf(xx, "%5.3f %8.4f %8.4f", para->hbLamda1, minEne, minRms);
            out << string(xx) << endl;
        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "hbLamda2") {
        for(double x= 0.1; x < 1.0;x = x*1.3) {
           
           para->hbLamda2 = x;
            et->updateAtomic(para);

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib,eiLib, et, inputFile, results, i-startID);
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
                results[0]->mergeResult(results[i]);
            }

            double minEne = results[0]->ene;
            double minRms = results[0]->rms;
            sprintf(xx, "%5.3f %8.4f %8.4f", x, minEne, minRms);
            out << string(xx) << endl;
        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "wtHb") {
        for(double x= 0.3; x < 5.0;x = x*1.3) {
           
           para->wtHb = x;
            et->updateAtomic(para);

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            cout << "start mp" << endl;
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, eiLib, et, inputFile,results, i-startID);
                jid++;
                thrPool->addTask(request);
            }
            while(true) {
                sleep(1);
                if(thrPool->getTaskCount() == 0) {
                    break;
                }
             }

            cout << "merge" << endl;
            for(int i=1;i<mp;i++){
                results[0]->mergeResult(results[i]);
            }

            
            double minEne = results[0]->ene;
            double minRms = results[0]->rms;

            cout << "print " << endl;
            sprintf(xx, "%5.3f %8.4f %8.4f", x, minEne, minRms);
            out << string(xx) << endl;
        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "wtBp1") {
        for(double x= 0.3; x < 6.0;x = x*1.3) {
           
           para->wtBp1 = x;
           et->bpET->updateForceFieldPara(para);

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, eiLib, et, inputFile, results, i-startID);
                jid++;
                thrPool->addTask(request);
            }
            while(true) {
                sleep(1);
                if(thrPool->getTaskCount() == 0) {
                    break;
                }
             }

             cout << "merge" << endl;
            for(int i=1;i<mp;i++){
                results[0]->mergeResult(results[i]);
            }

            double minEne = results[0]->ene;
            double minRms = results[0]->rms;
            sprintf(xx, "%5.3f %8.4f %8.4f", x, minEne, minRms);
            out << string(xx) << endl;
        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "wtBp2") {
        for(double x= 0.3; x < 6.0;x = x*1.3) {
           
           para->wtBp2 = x;
           et->bpET->updateForceFieldPara(para);

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 
  
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, eiLib, et, inputFile, results, i-startID);
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
                results[0]->mergeResult(results[i]);
            }

            double minEne = results[0]->ene;
            double minRms = results[0]->rms;
            sprintf(xx, "%5.3f %8.4f %8.4f", x, minEne, minRms);
            out << string(xx) << endl;
        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
        else if(tag == "wtBp") {
        for(double x= 0.3; x < 6.0;x = x*1.3) {
           para->wtBp1 = x;
           para->wtBp2 = x;
           et->bpET->updateForceFieldPara(para);

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 
  
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, eiLib, et, inputFile, results, i-startID);
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
                results[0]->mergeResult(results[i]);
            }

            double minEne = results[0]->ene;
            double minRms = results[0]->rms;
            sprintf(xx, "%5.3f %8.4f %8.4f", x, minEne, minRms);
            out << string(xx) << endl;
        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "wtO4O2C2Nb") {
        for(double x= 0.0; x < 2.0;x = x+0.2) {
           
           para->wtO4O2C2Nb = x;
            et->updateAtomic(para);

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 
 
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, eiLib, et, inputFile,results, i-startID);
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
                results[0]->mergeResult(results[i]);
            }

            double minEne = results[0]->ene;
            double minRms = results[0]->rms;
            sprintf(xx, "%5.3f %8.4f %8.4f", x, minEne, minRms);
            out << string(xx) << endl;
        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "wtO4O2C2Nnb") {
        for(double x= 0.0; x < 1.21;x = x+0.2) {
           
           para->wtO4O2C2Nnb = x;
            et->updateAtomic(para);

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 
   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, eiLib, et,  inputFile, results, i-startID);
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
                results[0]->mergeResult(results[i]);
            }

            double minEne = results[0]->ene;
            double minRms = results[0]->rms;
            sprintf(xx, "%5.3f %8.4f %8.4f", x, minEne, minRms);
            out << string(xx) << endl;
        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "wtBB") {
        //for(double x= 0.1; x < 2.0;x = x+0.1) {


        {   
            para->wtRibose = y;
            para->wtPho = y*0.7;
            et->updateAtomic(para);

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, eiLib, et, inputFile, results, i-startID);
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
                results[0]->mergeResult(results[i]);
            }

            double minEne = results[0]->ene;
            double minRms = results[0]->rms;
            sprintf(xx, "%5.3f %8.4f %8.4f", y, minEne, minRms);
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
 
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, eiLib, et, inputFile, results,  i-startID);
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
                results[0]->mergeResult(results[i]);
            }

            double minEne = results[0]->ene;
            double minRms = results[0]->rms;
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

            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testRefinement, moveLib, eiLib, et, inputFile, results, i-startID);
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
                results[0]->mergeResult(results[i]);
            }

            double minEne = results[0]->ene;
            double minRms = results[0]->rms;
            sprintf(xx, "%5.3f %8.4f %8.4f", x, minEne, minRms);
            out << string(xx) << endl;
        }
        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;        
    }

    else if(tag == "test") {
        int stepNum = testSamplingStepNum(moveLib, eiLib, et, inputFile);
        out << stepNum << endl;
    }

    out.close();
    delete moveLib;
    delete et;
    delete para;

}
