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



class singlePredictionResult{
public:
    int posNum;
    vector<int> posList;
    vector<double> rmsList;
    vector<double> eneList;

    singlePredictionResult(){
        this->posNum = 0;
    }

    double getMeanRMS(){
        double rms = 0;
        for(int i=0;i<rmsList.size();i++){
            rms += rmsList[i];
        }
        return rms/rmsList.size();
    }

    double getMeanEne() {
        double ene = 0;
        for(int i=0;i<eneList.size();i++){
            ene += eneList[i];
        }
        return ene/eneList.size();
    }

    void addPos(int pos, double ene, double rms){
        posList.push_back(pos);
        if(rms > 2.5)
            rms = 2.5;
        rmsList.push_back(rms);
        eneList.push_back(ene);
        this->posNum = posList.size();
    }

    void mergeResult(singlePredictionResult* other){
        if(this->posNum != other->posNum){
            cout << "pos num not consistent" << endl;
        }
        for(int i=0;i<posNum;i++){
            if(other->eneList[i] < this->eneList[i]){
                this->rmsList[i] = other->rmsList[i];
                this->eneList[i] = other->eneList[i];
            }
        }
    }

    void clear(){
        this->posNum = 0;
        this->posList.clear();
        this->rmsList.clear();
        this->eneList.clear();
    }
};


int testSingleBasePredictionPosList(NuPairMoveSetLibrary* moveLib, EdgeInformationLib* eiLib, RnaEnergyTable* et, const string& libType,  const string& inputFileList, singlePredictionResult** outList, int mpID){
	
    BasePairLib* pairLib = new BasePairLib(libType);
    RotamerLib* rotLib = new RotamerLib();
    AtomLib* atLib = new AtomLib();

    ifstream file;
    file.open(inputFileList.c_str(), ios::in);
    string line;

    vector<string> inputList;
    while(getline(file, line)){
        inputList.push_back(line);
    }

    for(int i=0;i<inputList.size();i++){
        NSPtools::InputParser input(inputList[i]);
        string cst = input.getValue("cst");
        int pos = -1;
        for(int k=0;k<cst.length();k++){
            if(cst[k] == '0')
                pos = k;
        }

        NuGraph* graph = new NuGraph(inputList[i], rotLib, atLib, pairLib, moveLib, eiLib, et);

        cout << "init for single residue prediction" << endl;
        graph->initForSingleResiduePrediction(inputList[i], pos);
	    graph->initRandWeight();
        
        cout << "new tree" << endl;
        NuTree* tree = new NuTree(graph);
        cout << "mst" << endl;
	    graph->MST_kruskal(tree);
	    tree->printEdges();
	    tree->updateNodeInfo(1.0, 1.0);
	    tree->updateEdgeInfo(1.0, 1.0);
	    tree->updateSamplingInfo();

	    graphInfo* gi = tree->runAtomicMC();
        double rms = gi->rmsd(graph->initInfo, pos);
        gi->setRMS(rms);

        outList[mpID]->addPos(pos, gi->ene, gi->rms);
       
        delete gi;
	    delete tree;
	    delete graph;
    }

    delete pairLib;
    delete rotLib;
    delete atLib;

    return 0;
}

int main(int argc, char** argv){
        
    if(argc == 1) {
        cout << "trainSingle2 $inputList $list -lib $LIBTYPE -tag $tag" << endl;
        return EXIT_SUCCESS;
    }
    CmdArgs cmdArgs{argc, argv};


    string inputFileList = cmdArgs.getValue("-list");
    string outputFile = cmdArgs.getValue("-out");

    string tag = cmdArgs.getValue("-tag");
    string libType = "stat";
    int pairType = atoi(cmdArgs.getValue("-p").c_str());
    int clusterType = atoi(cmdArgs.getValue("-c").c_str());


    cout << "load move lib: " << endl;
	NuPairMoveSetLibrary* moveLib = new NuPairMoveSetLibrary(libType, true, 1);
	moveLib->load();

    EdgeInformationLib* eiLib = new EdgeInformationLib();


    string outpdb;

    ofstream out;


    char xx[200];
    int mp = 8;
    int startID = 0;
	clock_t start = clock();

    if(cmdArgs.specifiedOption("-mp"))
        mp = atoi(cmdArgs.getValue("-mp").c_str());

    singlePredictionResult** outList = new singlePredictionResult*[mp];
    for(int i=0;i<mp;i++){
        outList[i] = new singlePredictionResult();
    }

    ForceFieldPara* para = new ForceFieldPara();

    para->kStepNum1 = 200;
    para->kStepNum2 = 300;
    para->kStepNum3 = 600;
    para->kNodeFreq = 0.8;

    para->libType = libType;

    srand(time(0));

    cout << "load energy table" << endl;
    RnaEnergyTable* et = new RnaEnergyTable(para);
	et->loadAtomicEnergy();


    cout << "task: " << tag << endl;
    cout << "a" << endl;

    if(tag == "clashNb") {
        out.open(outputFile.c_str(), ios::out);
        for(double x= 1.0; x < 10.0;x = x*1.3) {
           
           para->clashLamdaNb = x;
           et->updateAtomic(para);

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outList[i]->clear();
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePredictionPosList, moveLib, eiLib, et, libType, inputFileList, outList, i-startID);
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
                outList[0]->mergeResult(outList[i]);
            }

            double minEne = 0.0;
            double minRMS = 0.0;
            for(int i=0;i<outList[0]->posNum;i++){
                minEne += outList[0]->eneList[i]/outList[0]->posNum;
                minRMS += outList[0]->rmsList[i]/outList[0]->posNum;
            }
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRMS, outList[0]->posNum);
            out << string(xx) << endl;

        }

        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "clashNnb") {
        out.open(outputFile.c_str(), ios::out);
        for(double x= 1.0; x < 10.0;x = x*1.3) {
           
            para->clashLamdaNnb = x;
            et->updateAtomic(para);

             shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outList[i]->clear();
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                 request->asynBind(testSingleBasePredictionPosList, moveLib, eiLib, et, libType, inputFileList, outList, i-startID);
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
                outList[0]->mergeResult(outList[i]);
            }

            double minEne = 0.0;
            double minRMS = 0.0;
            for(int i=0;i<outList[0]->posNum;i++){
                minEne += outList[0]->eneList[i]/outList[0]->posNum;
                minRMS += outList[0]->rmsList[i]/outList[0]->posNum;
            }
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRMS, outList[0]->posNum);
            out << string(xx) << endl;

        }

        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "hbLamda1") {
        out.open(outputFile.c_str(), ios::out);
        for(double x= 0.03; x < 0.5;x = x*1.3) {
           
            para->hbLamda1 = x;
            et->updateAtomic(para);

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outList[i]->clear();
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePredictionPosList, moveLib, eiLib, et, libType, inputFileList, outList, i-startID);
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
                outList[0]->mergeResult(outList[i]);
            }

            double minEne = 0.0;
            double minRMS = 0.0;
            for(int i=0;i<outList[0]->posNum;i++){
                minEne += outList[0]->eneList[i]/outList[0]->posNum;
                minRMS += outList[0]->rmsList[i]/outList[0]->posNum;
            }
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRMS, outList[0]->posNum);
            out << string(xx) << endl;

        }

        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "hbLamda2") {
        out.open(outputFile.c_str(), ios::out);
        for(double x= 0.1; x < 1.0;x = x*1.3) {
           
            para->hbLamda2 = x;
            et->updateAtomic(para);

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outList[i]->clear();
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePredictionPosList, moveLib, eiLib, et, libType, inputFileList, outList, i-startID);
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
                outList[0]->mergeResult(outList[i]);
            }

            double minEne = 0.0;
            double minRMS = 0.0;
            for(int i=0;i<outList[0]->posNum;i++){
                minEne += outList[0]->eneList[i]/outList[0]->posNum;
                minRMS += outList[0]->rmsList[i]/outList[0]->posNum;
            }
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRMS, outList[0]->posNum);
            out << string(xx) << endl;

        }

        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "wtHb") {
        out.open(outputFile.c_str(), ios::out);
        for(double x= 0.3; x < 5.0;x = x*1.3) {
           
            para->wtHb = x;
            et->updateAtomic(para);

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outList[i]->clear();
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
               request->asynBind(testSingleBasePredictionPosList, moveLib, eiLib, et, libType, inputFileList, outList, i-startID);
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
                outList[0]->mergeResult(outList[i]);
            }

            double minEne = 0.0;
            double minRMS = 0.0;
            for(int i=0;i<outList[0]->posNum;i++){
                minEne += outList[0]->eneList[i]/outList[0]->posNum;
                minRMS += outList[0]->rmsList[i]/outList[0]->posNum;
            }
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRMS, outList[0]->posNum);
            out << string(xx) << endl;

        }

        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "wtBp1") {
        out.open(outputFile.c_str(), ios::out);
        for(double x= 0.3; x < 6.0;x = x*1.3) {
           
           para->wtBp1 = x;
           et->bpET->updateForceFieldPara(para);

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outList[i]->clear();
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                 request->asynBind(testSingleBasePredictionPosList, moveLib, eiLib, et, libType, inputFileList, outList, i-startID);
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
                outList[0]->mergeResult(outList[i]);
            }

            double minEne = 0.0;
            double minRMS = 0.0;
            for(int i=0;i<outList[0]->posNum;i++){
                minEne += outList[0]->eneList[i]/outList[0]->posNum;
                minRMS += outList[0]->rmsList[i]/outList[0]->posNum;
            }
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRMS, outList[0]->posNum);
            out << string(xx) << endl;

        }

        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "wtBp2") {
        out.open(outputFile.c_str(), ios::out);
        for(double x= 0.3; x < 6.0;x = x*1.3) {
           
           para->wtBp2 = x;
            et->bpET->updateForceFieldPara(para);
            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outList[i]->clear();
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
               request->asynBind(testSingleBasePredictionPosList, moveLib, eiLib, et, libType, inputFileList, outList, i-startID);
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
                outList[0]->mergeResult(outList[i]);
            }

            double minEne = 0.0;
            double minRMS = 0.0;
            for(int i=0;i<outList[0]->posNum;i++){
                minEne += outList[0]->eneList[i]/outList[0]->posNum;
                minRMS += outList[0]->rmsList[i]/outList[0]->posNum;
            }
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRMS, outList[0]->posNum);
            out << string(xx) << endl;

        }

        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "wtO4O2C2Nb") {
        out.open(outputFile.c_str(), ios::out);
        for(double x= 0.0; x < 2.0;x = x+0.2) {
           
            para->wtO4O2C2Nb = x;
            et->updateAtomic(para);

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outList[i]->clear();
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePredictionPosList, moveLib, eiLib, et, libType, inputFileList, outList, i-startID);
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
                outList[0]->mergeResult(outList[i]);
            }

            double minEne = 0.0;
            double minRMS = 0.0;
            for(int i=0;i<outList[0]->posNum;i++){
                minEne += outList[0]->eneList[i]/outList[0]->posNum;
                minRMS += outList[0]->rmsList[i]/outList[0]->posNum;
            }
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRMS, outList[0]->posNum);
            out << string(xx) << endl;

        }

        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "wtO4O2C2Nnb") {
        out.open(outputFile.c_str(), ios::out);
        for(double x= 0.0; x < 1.21;x = x+0.2) {
           
            para->wtO4O2C2Nnb = x;
            et->updateAtomic(para);

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outList[i]->clear();
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
               request->asynBind(testSingleBasePredictionPosList, moveLib, eiLib, et, libType, inputFileList, outList, i-startID);
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
                outList[0]->mergeResult(outList[i]);
            }

            double minEne = 0.0;
            double minRMS = 0.0;
            for(int i=0;i<outList[0]->posNum;i++){
                minEne += outList[0]->eneList[i]/outList[0]->posNum;
                minRMS += outList[0]->rmsList[i]/outList[0]->posNum;
            }
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRMS, outList[0]->posNum);
            out << string(xx) << endl;

        }

        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "wtOxy") {
        out.open(outputFile.c_str(), ios::out);
        for(double x= 0.1; x < 4.0;x = x+0.2) {
            for(int k=0;k<16;k++) {
                para->wtRiboseOxy[k] = x;
            }
            
            delete et->roET;
            et->roET = new RiboseOxygenEnergyTable(para);

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outList[i]->clear();
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePredictionPosList, moveLib, eiLib, et, libType, inputFileList, outList, i-startID);
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
                outList[0]->mergeResult(outList[i]);
            }

            double minEne = 0.0;
            double minRMS = 0.0;
            for(int i=0;i<outList[0]->posNum;i++){
                minEne += outList[0]->eneList[i]/outList[0]->posNum;
                minRMS += outList[0]->rmsList[i]/outList[0]->posNum;
            }
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRMS, outList[0]->posNum);
            out << string(xx) << endl;

        }

        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;
    }
    else if(tag == "T0"){
        out.open(outputFile.c_str(), ios::out);
        for(double x= 2.0; x < 50.0;x = x*1.3) {
            para->T0 = x;
            para->T1 = para->T0*0.1;
            para->T2 = para->T1*0.1;
            para->T3 = para->T2*0.1;

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outList[i]->clear();
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePredictionPosList, moveLib, eiLib, et, libType, inputFileList, outList, i-startID);
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
                outList[0]->mergeResult(outList[i]);
            }

            double minEne = 0.0;
            double minRMS = 0.0;
            for(int i=0;i<outList[0]->posNum;i++){
                minEne += outList[0]->eneList[i]/outList[0]->posNum;
                minRMS += outList[0]->rmsList[i]/outList[0]->posNum;
            }
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRMS, outList[0]->posNum);
            out << string(xx) << endl;

        }

        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;       
    }
    else if(tag == "step1"){
        out.open(outputFile.c_str(), ios::out);
    
        for(double x= 50; x < 1000;x = x*1.5) {
           
            para->kStepNum1 = (int)x;
            para->kStepNum2 = (int)((1000-x)*0.5);
            para->kStepNum3 = (int)((1000-x)*0.5);

            cout << "stepNum: " <<  para->kStepNum1 << " " << para->kStepNum2 << " " << para->kStepNum3 << endl;

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outList[i]->clear();
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
               request->asynBind(testSingleBasePredictionPosList, moveLib, eiLib, et, libType, inputFileList, outList, i-startID);
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
                outList[0]->mergeResult(outList[i]);
            }

            double minEne = 0.0;
            double minRMS = 0.0;
            for(int i=0;i<outList[0]->posNum;i++){
                minEne += outList[0]->eneList[i]/outList[0]->posNum;
                minRMS += outList[0]->rmsList[i]/outList[0]->posNum;
            }
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRMS, outList[0]->posNum);
            out << string(xx) << endl;

        }
    }
    else if(tag == "step2"){
        out.open(outputFile.c_str(), ios::out);
    
        for(double x= 50; x < 1000;x = x*1.5) {
           
            para->kStepNum2 = (int)x;
            para->kStepNum1 = (int)((1000-x)*0.5);
            para->kStepNum3 = (int)((1000-x)*0.5);

            cout << "stepNum: " <<  para->kStepNum1 << " " << para->kStepNum2 << " " << para->kStepNum3 << endl;


            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outList[i]->clear();
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePredictionPosList, moveLib, eiLib, et, libType, inputFileList, outList, i-startID);
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
                outList[0]->mergeResult(outList[i]);
            }

            double minEne = 0.0;
            double minRMS = 0.0;
            for(int i=0;i<outList[0]->posNum;i++){
                minEne += outList[0]->eneList[i]/outList[0]->posNum;
                minRMS += outList[0]->rmsList[i]/outList[0]->posNum;
            }
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRMS, outList[0]->posNum);
            out << string(xx) << endl;

        }

        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;       
    }
    else if(tag == "step3"){
        out.open(outputFile.c_str(), ios::out);
    
        for(double x= 50; x < 1000;x = x*1.5) {
           
            para->kStepNum3 = (int)x;
            para->kStepNum1 = (int)((1000-x)*0.5);
            para->kStepNum2 = (int)((1000-x)*0.5);

            cout << "stepNum: " <<  para->kStepNum1 << " " << para->kStepNum2 << " " << para->kStepNum3 << endl;

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outList[i]->clear();
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
               request->asynBind(testSingleBasePredictionPosList, moveLib, eiLib, et, libType, inputFileList, outList, i-startID);
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
                outList[0]->mergeResult(outList[i]);
            }

            double minEne = 0.0;
            double minRMS = 0.0;
            for(int i=0;i<outList[0]->posNum;i++){
                minEne += outList[0]->eneList[i]/outList[0]->posNum;
                minRMS += outList[0]->rmsList[i]/outList[0]->posNum;
            }
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRMS, outList[0]->posNum);
            out << string(xx) << endl;

        }
   }
    else if(tag == "kNode"){
        out.open(outputFile.c_str(), ios::out);
    
        for(double x= 0.1; x < 5.0;x = x*1.5) {
           
            para->kStepNum3 = 333;
            para->kStepNum1 = 333;
            para->kStepNum2 = 333;

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outList[i]->clear();
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePredictionPosList, moveLib, eiLib, et, libType, inputFileList, outList, i-startID);
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
                outList[0]->mergeResult(outList[i]);
            }

            double minEne = 0.0;
            double minRMS = 0.0;
            for(int i=0;i<outList[0]->posNum;i++){
                minEne += outList[0]->eneList[i]/outList[0]->posNum;
                minRMS += outList[0]->rmsList[i]/outList[0]->posNum;
            }
            
            sprintf(xx, "%4.2f %8.4f %8.4f %d", x, minEne, minRMS, outList[0]->posNum);
            out << string(xx) << endl;

        }
    }
    else if(tag == "bw") {
        out.open(outputFile.c_str(), ios::out);
        {
            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outList[i]->clear();
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
               request->asynBind(testSingleBasePredictionPosList, moveLib, eiLib, et, libType, inputFileList, outList, i-startID);
                jid++;
                thrPool->addTask(request);
            }

            while(true) {
                sleep(1);
                if(thrPool->getTaskCount() == 0) {
                    cout << "finish mc" << endl;
                    break;
                }
             }

            cout << "merge result" << endl;


            for(int i=1;i<mp;i++){
                cout << "mp: " << i  << endl;
                outList[0]->mergeResult(outList[i]);
            }

            cout << "posNum: " << outList[0]->posNum << endl;

            double minEne = 0.0;
            double minRMS = 0.0;
            for(int i=0;i<outList[0]->posNum;i++){
                minEne += outList[0]->eneList[i]/outList[0]->posNum;
                minRMS += outList[0]->rmsList[i]/outList[0]->posNum;
            }

            cout << minEne << " " << minRMS << endl;
            
            sprintf(xx, "%4.2f %8.4f %d", minEne, minRMS, outList[0]->posNum);
            out << string(xx) << endl;
        }
    }
    else if(tag == "test") {
        for(int x=0;x<3;x++) {
            ofstream out2;
            sprintf(xx, "%s-%d", outputFile.c_str(), x);
            string fileName2 = string(xx);
            out2.open(fileName2.c_str(), ios::out);

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outList[i]->clear();
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePredictionPosList, moveLib, eiLib, et, libType, inputFileList, outList, i-startID);
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
                outList[0]->mergeResult(outList[i]);
            }

            double minEne = 0.0;
            double minRMS = 0.0;
            for(int i=0;i<outList[0]->posNum;i++){
                sprintf(xx, "%d %8.4f %8.4f", i, outList[0]->rmsList[i], outList[0]->eneList[i]);
                out2 << string(xx) << endl;
            }
            out2.close();
        }

        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;        
    }
    else if(tag == "subEne") {

        out.open(outputFile.c_str(), ios::out);

        for(int x=0;x<10;x++) {
            et->bpET->addSubEnergy(pairType, clusterType, x*0.1);

           

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outList[i]->clear();
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePredictionPosList, moveLib, eiLib, et, libType, inputFileList, outList, i-startID);
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
                outList[0]->mergeResult(outList[i]);
            }

            double meanEne = outList[0]->getMeanEne();
            double meanRms = outList[0]->getMeanRMS();
            out << x*0.1 << " " << meanRms << " " << meanEne << endl;
     
        }
        out.close();

        clock_t end1 = clock();
	    cout << "mp: " << mp <<" " << "time: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;        
    }

    delete moveLib;
    delete et;
    delete para;

    for(int i=0;i<mp;i++){
        delete outList[i];
    }
    delete [] outList;

    

}
