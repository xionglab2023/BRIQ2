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

int testSamplingStepNum(NuPairMoveSetLibrary* moveLib, EdgeMoveClustersLib* eiLib, RnaEnergyTable* et, const string& inputFile){
    BasePairLib* pairLib = new BasePairLib();
	RotamerLib* rotLib = new RotamerLib();
	AtomLib* atLib = new AtomLib();
	NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, eiLib, et);
	graph->generateRandomEdgePartition(2);
	SamplingGraph* tree = new SamplingGraph(graph);
    tree->updatePartitionInfo();
    tree->updateSamplingInfo();  
    return (int)(tree->totalSamp);
}

int testSingleBasePrediction(NuPairMoveSetLibrary* moveLib, EdgeMoveClustersLib* eiLib, RnaEnergyTable* et, BasePairLib* pairLib,RotamerLib* rotLib,AtomLib* atLib, const string& inputFile, singlePredictionResult** outList, int mpID){
	
    NSPtools::InputParser input(inputFile);
    string cst = input.getValue("cst");
    char xx[200];

    for(int i=0;i<cst.length();i++){
        
        if(cst[i] == 'F') continue;
        sprintf(xx, "%s-%d", inputFile.c_str(), i);
        string singlePosInput = string(xx);
    
        cout << "pos: " << i << " " << cst[i] << endl;
        cout << "init graph" << endl;

        NSPtools::InputParser input2(singlePosInput);
        string cst2 = input2.getValue("cst");
        int pos = -1;
        for(int k=0;k<cst2.length();k++){
            if(cst2[k] == '0')
                pos = k;
        }

        NuGraph* graph = new NuGraph(singlePosInput, rotLib, atLib, pairLib, moveLib, eiLib, et);

        cout << "init for single residue prediction" << endl;
        graph->initForSingleResiduePrediction(singlePosInput, pos);

        graph->generateRandomEdgePartition(2);
        SamplingGraph* tree = new SamplingGraph(graph);
        tree->updatePartitionInfo();
        tree->updateSamplingInfo();
	    
        cout << "run mc" << endl;
	    graphInfo* gi = tree->runAtomicMC();
        double rms = gi->rmsd(graph->initInfo, pos);
        gi->setRMS(rms);

        outList[mpID]->addPos(pos, gi->ene, gi->rms);
       
        delete gi;
	    delete tree;
	    delete graph;
    }


    return 0;
}

int testSingleBasePredictionPosList(NuPairMoveSetLibrary* moveLib, EdgeMoveClustersLib* eiLib, RnaEnergyTable* et, BasePairLib* pairLib,RotamerLib* rotLib,AtomLib* atLib, const string& inputFileList, singlePredictionResult** outList, int mpID){
	

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
        graph->generateRandomEdgePartition(2);
        SamplingGraph* tree = new SamplingGraph(graph);
        tree->updatePartitionInfo();
        tree->updateSamplingInfo();

	    graphInfo* gi = tree->runAtomicMC();
        double rms = gi->rmsd(graph->initInfo, pos);
        gi->setRMS(rms);

        outList[mpID]->addPos(pos, gi->ene, gi->rms);
       
        delete gi;
	    delete tree;
	    delete graph;
    }

    return 0;
}

int testSingleBasePredictionPrintPDB(NuPairMoveSetLibrary* moveLib,EdgeMoveClustersLib* eiLib, RnaEnergyTable* et, BasePairLib* pairLib,RotamerLib* rotLib,AtomLib* atLib, const string& inputFile, ofstream& out, const string& pdbFile, const string& pdbID){
	
    NSPtools::InputParser input(inputFile);
    string cst = input.getValue("cst");


    string pdbTag = pdbFile.substr(0, pdbFile.length()-4);
    char xx[200];
    for(int i=0;i<cst.length();i++){
        if(cst[i] == 'F') continue;

        sprintf(xx, "%s-%d.pdb", pdbTag.c_str(), i);
        string pdbout = string(xx);

        sprintf(xx, "%s-%d", inputFile.c_str(), i);
        string singlePosInput = string(xx);

        NSPtools::InputParser input2(singlePosInput);
        string cst2 = input2.getValue("cst");
        int pos = -1;
        for(int k=0;k<cst2.length();k++){
            if(cst2[k] == '0')
                pos = k;
        }

        NuGraph* graph = new NuGraph(singlePosInput, rotLib, atLib, pairLib, moveLib, eiLib, et);
        graph->initForSingleResiduePrediction(singlePosInput, pos);
        graph->generateRandomEdgePartition(2);
        SamplingGraph* tree = new SamplingGraph(graph);
        tree->updatePartitionInfo();
        tree->updateSamplingInfo();
        
	    graphInfo* gi = tree->runAtomicMC();
        double rms = gi->rmsd(graph->initInfo, pos);
        out  << pdbID << " " <<  i << " " <<  gi->ene  << " " << rms << endl;
        gi->printPDB(pdbout);
        //graph->printEnergy();

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

    string libType = "stat";

    cout << "load move lib: " << endl;
	NuPairMoveSetLibrary* moveLib = new NuPairMoveSetLibrary(libType, true, 1);
	moveLib->load();

    BasePairLib* pairLib = new BasePairLib(libType);
    EdgeMoveClustersLib* eiLib = new EdgeMoveClustersLib();

    string inputFile = cmdArgs.getValue("-in");
    string tag = cmdArgs.getValue("-tag");
    string outputFile = cmdArgs.getValue("-out");
    //string bw = cmdArgs.getValue("-bw");
    string outpdb;
    if(cmdArgs.specifiedOption("-outpdb"))
        outpdb = cmdArgs.getValue("-outpdb");
    string pdbID = "pdb";
     if(cmdArgs.specifiedOption("-pdbID"))
        pdbID = cmdArgs.getValue("-pdbID");

    ofstream out;
    out.open(outputFile.c_str(), ios::out);

    char xx[200];
    int mp = 16;
    int startID = 0;
	clock_t start = clock();

    singlePredictionResult** outList = new singlePredictionResult*[mp];
    for(int i=0;i<mp;i++){
        outList[i] = new singlePredictionResult();
    }

    ForceFieldPara* para = new ForceFieldPara();
   
    para->libType = libType;

    para->kStepNum1 = 100;
    para->kStepNum2 = 300;
    para->kStepNum3 = 600;
    para->kNodeFreq = 0.8;

    srand(time(0));

    cout << "load energy table" << endl;
    RnaEnergyTable* et = new RnaEnergyTable(para);
	et->loadAtomicEnergy();

    cout << "init pairLib rotLib atLib" << endl;
   
	RotamerLib* rotLib = new RotamerLib();
	AtomLib* atLib = new AtomLib();


    cout << "task: " << tag << endl;
    cout << "a" << endl;

    if(tag == "clashNb") {
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
                request->asynBind(testSingleBasePrediction, moveLib, eiLib, et, pairLib, rotLib, atLib, inputFile, outList, i-startID);
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
                request->asynBind(testSingleBasePrediction, moveLib,  eiLib, et, pairLib, rotLib, atLib, inputFile, outList, i-startID);
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
                request->asynBind(testSingleBasePrediction, moveLib, eiLib, et, pairLib, rotLib, atLib, inputFile, outList, i-startID);
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
                request->asynBind(testSingleBasePrediction, moveLib, eiLib, et, pairLib, rotLib, atLib, inputFile, outList, i-startID);
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
                request->asynBind(testSingleBasePrediction, moveLib, eiLib, et, pairLib, rotLib, atLib, inputFile, outList, i-startID);
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
    else if(tag == "wtBp") {
        for(double x= 0.5; x < 5.0;x = x*1.15) {
           
           para->wtBp1 = x;
           para->wtBp2 = x;
           et->bpET->updateForceFieldPara(para);

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outList[i]->clear();
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePrediction, moveLib, eiLib, et, pairLib, rotLib, atLib, inputFile, outList, i-startID);
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
                request->asynBind(testSingleBasePrediction, moveLib, eiLib, et, pairLib, rotLib, atLib, inputFile, outList, i-startID);
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
                request->asynBind(testSingleBasePrediction, moveLib, eiLib, et, pairLib, rotLib, atLib, inputFile, outList, i-startID);
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
                request->asynBind(testSingleBasePrediction, moveLib, eiLib, et, pairLib, rotLib, atLib, inputFile, outList, i-startID);
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
                request->asynBind(testSingleBasePrediction, moveLib, eiLib, et, pairLib, rotLib, atLib, inputFile, outList, i-startID);
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
    else if(tag == "wtBB") {
        
        for(double x= 0.1; x < 2.0;x = x+0.1) {
           
            para->wtRibose = x;
            para->wtPho = x*0.7;
            et->updateAtomic(para);

            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outList[i]->clear();
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePrediction, moveLib,eiLib, et, pairLib, rotLib, atLib, inputFile, outList, i-startID);
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
                request->asynBind(testSingleBasePrediction, moveLib, eiLib, et, pairLib, rotLib, atLib, inputFile, outList, i-startID);
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
                request->asynBind(testSingleBasePrediction, moveLib, eiLib, et, pairLib, rotLib, atLib, inputFile, outList, i-startID);
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
                request->asynBind(testSingleBasePrediction, moveLib, eiLib, et, pairLib, rotLib, atLib, inputFile, outList, i-startID);
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
                request->asynBind(testSingleBasePrediction, moveLib, eiLib, et, pairLib, rotLib, atLib, inputFile, outList, i-startID);
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
                request->asynBind(testSingleBasePrediction, moveLib, eiLib, et, pairLib, rotLib, atLib, inputFile, outList, i-startID);
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
                request->asynBind(testSingleBasePrediction, moveLib, eiLib, et, pairLib, rotLib, atLib, inputFile, outList, i-startID);
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
        {
         
            shared_ptr<ThreadPool> thrPool(new ThreadPool(mp));
            size_t jid = 0; 

            for(int i=0;i<mp;i++){
                outList[i]->clear();
            }   
            for(int i=startID;i<startID+mp;i++) {
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(testSingleBasePrediction, moveLib, eiLib, et, pairLib, rotLib, atLib, inputFile, outList, i-startID);
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
    else if(tag == "pdb") {
        testSingleBasePredictionPrintPDB(moveLib, eiLib, et, pairLib, rotLib, atLib, inputFile, out, outpdb, pdbID);
    }

    out.close();
    delete moveLib;
    delete pairLib;
    delete eiLib;
    delete rotLib;
    delete atLib;
    delete et;
    delete para;

    for(int i=0;i<mp;i++){
        delete outList[i];
    }
    delete [] outList;

}
