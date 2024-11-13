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

int runRefinementsub(NuPairMoveSetLibrary* moveLib, EdgeInformationLib* eiLib, RnaEnergyTable* et, const string& inputFile, const string& outFile, int randSeed, double kStep, bool printDetail){

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

	int* subGraphPosList = new int[graph->seqLen]();
	int* fixedPositions = new int[graph->seqLen]();
	vector<int> outsubGraphPosList;
	vector<int> outupdatePosList;
	vector<vector<int>> Bclusters;
	vector<vector<int>> Cclusters;
	int corePos = 0;
	int NumberOfCycles = 0;//记录循环次数
	vector<int> count(graph->seqLen, 0);  // 初始化每个元素的计数
	NuGraph* subGraph = nullptr;

	while (true)
	{	
		corePos = rand() % (graph->seqLen);  // 随机选择一个位置作为 corePos

		subGraph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, eiLib, et);
		graph->generateSubGraph(inputFile, corePos, subGraphPosList, fixedPositions, subGraph, outsubGraphPosList, outupdatePosList, Bclusters, Cclusters);

		//在run MC之前/更新节点坐标前记录cm关系
		vector<vector<CsMove>> csMoveLists;
		for(size_t i = 0; i < Bclusters.size(); i++){
			vector<CsMove> csMoveList;
			int Bpos = Bclusters[i][0];
			// cout << "Bpos:" << Bpos << endl;
			for(size_t j = 0; j < Cclusters[i].size(); j++) {
				int Cpos = Cclusters[i][j];
				// if (j ==0) {cout << "Cpos:" << endl;}
				csMoveList.push_back(graph->allNodes[Cpos]->baseConf->cs1 - graph->allNodes[Bpos]->baseConf->cs1);
				// cout << j << "-" << Cpos << " ";
			}
			cout << endl;
			csMoveLists.push_back(csMoveList);
		}

		subGraph->initRandWeight();
		NuTree* subtree = new NuTree(subGraph);
		subGraph->MST_kruskal(subtree);
		subtree->updateNodeInfo(1.0, 1.0);
		subtree->updateEdgeInfo(1.0, 1.0);
		subtree->updateSamplingInfo();

		cout << "run MC" << endl;
		graphInfo* subgi = subtree->runAtomicMC();

		// const string& outputFile1 = "outFile_1105_3.pdb";
		// subgi->printPDB(outputFile1);

		cout << "update graph->allNodes" << endl;
		for (size_t i = 0; i < outsubGraphPosList.size(); i++) {
			LocalFrame cs1 = subGraph->allNodes[i]->baseConf->cs1;
			// RiboseRotamer* rot = new RiboseRotamer(*subGraph->allNodes[i]->riboseConf->rot);

			graph->allNodes[outsubGraphPosList[i]]->updateCoordinate(cs1);
			graph->allNodes[outsubGraphPosList[i]]->acceptCoordMove();

			// graph->allNodes[outsubGraphPosList[i]]->riboseConfTmp->updateRotamer(rot);
			graph->allNodes[outsubGraphPosList[i]]->riboseConfTmp->rot->copyValueFrom(subGraph->allNodes[i]->riboseConf->rot);
			graph->allNodes[outsubGraphPosList[i]]->riboseConfTmp->updateRotamer(graph->allNodes[outsubGraphPosList[i]]->riboseConfTmp->rot);
			graph->allNodes[outsubGraphPosList[i]]->riboseConf->copyValueFrom(graph->allNodes[outsubGraphPosList[i]]->riboseConfTmp);	
			// delete rot;
			// cout << "node " << outsubGraphPosList[i] << " updated" << endl;
		}	

		for(size_t i = 0; i < Bclusters.size(); i++){
			int Bpos = Bclusters[i][0];
			// cout << "Bpos:" << Bpos << endl;
			for(size_t j = 0; j < Cclusters[i].size(); j++) {
				int Cpos = Cclusters[i][j];
				// cout << "Cpos:" << Cpos << " ";
				LocalFrame cs1 = graph->allNodes[Bpos]->baseConfTmp->cs1 + csMoveLists[i][j];
				graph->allNodes[Cpos]->updateCoordinate(cs1);
				graph->allNodes[Cpos]->acceptCoordMove();
			}
			// cout << endl;
		}	
		graph->initPho();// 更新磷酸基团坐标

		vector<int> subconnectionBreakPoints;
		for(size_t i = 0; i < graph->seqLen; i++){
			if(graph->connectToDownstream[i]) { // 首先判断 connectToDownstream[i] 是否为 true
				// 检查条件：i 在 outsubGraphPosList 中且 i+1 不在，或者 i 不在且 i+1 在
				bool iInOutsubGraph = (find(outsubGraphPosList.begin(), outsubGraphPosList.end(), i) != outsubGraphPosList.end());
				bool iPlusOneInOutsubGraph = (find(outsubGraphPosList.begin(), outsubGraphPosList.end(), i + 1) != outsubGraphPosList.end());

				if ((iInOutsubGraph && !iPlusOneInOutsubGraph) || (!iInOutsubGraph && iPlusOneInOutsubGraph)) {
					subconnectionBreakPoints.push_back(i);
				}
			}
		}

		cout << "update phosphate: ";
		for(size_t i = 0;i < subconnectionBreakPoints.size();i++){
			int j = subconnectionBreakPoints[i];
			// cout << "subconnectionBreakPoints[" << i << "] = " <<  subconnectionBreakPoints[i] << ' ';
			graph->et->pb->buildPhosphate(graph->allNodes[j]->riboseConfTmp, graph->allNodes[j+1]->riboseConfTmp, graph->allNodes[j]->phoConfTmp);
			graph->allNodes[j]->phoConf->copyValueFrom(graph->allNodes[j]->phoConfTmp);
		}
		cout << endl;

		cout << "update graph->allEdges";
		for(size_t i = 0; i < graph->seqLen; i++){
			for(size_t j = 0; j < graph->seqLen; j++){
				graph->allEdges[i*graph->seqLen +j]->cm = graph->allNodes[j]->baseConf->cs1 - graph->allNodes[i]->baseConf->cs1;
				graph->allEdges[i*graph->seqLen +j]->cmTmp = graph->allEdges[i*graph->seqLen +j]->cm;
				// cout << j << "-" << i << " " ;
			}
		}
		cout << endl;

		vector<int> ().swap(subconnectionBreakPoints);
		vector<vector<CsMove>> ().swap(csMoveLists);
		vector<vector<int>> ().swap(Bclusters);
		vector<vector<int>> ().swap(Cclusters);
		vector<int> ().swap(outsubGraphPosList);

		// 更新 count 数组，记录 subGraphPosList 中的元素出现次数
        for (int pos : outupdatePosList) {
            count[pos]++;
        }
        vector<int> ().swap(outupdatePosList);  // 清空 outsubGraphPosList

        // 检查是否所有元素都出现至少三次
        bool allElementsMet = true;
        for (size_t i = 0; i < graph->seqLen; i++) {
            if (count[i] < 3) {
                allElementsMet = false;
                break;
            }
        }
	
		delete subgi;  // 循环内删除 subgi
		subgi = nullptr;
		delete subtree;
		subtree = nullptr;
		delete subGraph;// 保留该信息，给主图输出使用  可以再看
        subGraph = nullptr;

		// if (NumberOfCycles == 1){
		// 	break;
		// }
		
        // 如果所有元素都至少出现三次，结束循环
        if (allElementsMet) {
			// 输出已优化的节点坐标
			cout << "optimized node times: " << endl;
			for (size_t j = 0; j < graph->seqLen; j++) {
				cout << count[j] << " ";
			}
			cout << endl;
            break;
        }else{
            
        }

        cout << "Number Of Cycles:" << NumberOfCycles << endl;
        NumberOfCycles++;
    }
    delete [] subGraphPosList;
    delete [] fixedPositions;
    count.clear();
	
	graphInfo* gi = graph->getGraphInfo();	
	gi->printPDB(outFile);

	if(printDetail){
		cout << "detail energy: " << endl;
		graph->printEnergy();
	}

	delete gi;
    delete pairLib;
    delete rotLib;
    delete atLib;
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
		// runRefinementFromPDB(moveLib, eiLib, et, pdbFile, outputFile, seed, kStep, cnt);
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
			// runRefinement(moveLib, eiLib, et, inputFile, outputFile, seed, kStep, printEnergyDetail);
			runRefinementsub(moveLib, eiLib,et, inputFile, outputFile, seed, kStep, printEnergyDetail);
		}
		else if(task == "predict" && key.length() == 0) {
			et->loadAtomicEnergy();
			// runPredict(moveLib, eiLib, et, inputFile, outputFile, seed, kStep);
		}
		else if(task == "predict"){
			int modelNum = 10;
			if(cmdArgs.specifiedOption("-n"))
				modelNum = atoi(cmdArgs.getValue("-n").c_str());

			et->loadAtomicEnergy();
			// cgKeyToAllAtomPDB(moveLib, eiLib, et, inputFile, outputFile, modelNum, seed, kStep);
		}
		else if(task == "cgModeling"){

			int modelNum = 10;
			int mp = 1;
			if(cmdArgs.specifiedOption("-n"))
				modelNum = atoi(cmdArgs.getValue("-n").c_str());

			et->loadCoarseGrainedEnergy();
			// runCGMC(moveLib, et, eiLib, inputFile, outputFile, modelNum, seed, kStep);


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


