/*
 * NuGraph.h
 *
 *  Created on: 2023��11��15��
 *      Author: nuc
 */

#ifndef PREDNA_NUGRAPH_H_
#define PREDNA_NUGRAPH_H_
// #define DEBUG

#include <vector>
#include <string>
#include <map>
#include <queue>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include "model/StructureModel.h"
#include "model/AtomLib.h"
#include "model/RotamerLib.h"
#include "model/BaseRotamer.h"
#include "model/RiboseRotamer.h"
#include "model/PhosphateRotamer.h"
#include "model/BasePairLib.h"
#include "model/AssignRNASS.h"
#include "geometry/RMSD.h"
#include "predNA/NuMoveSet.h"
#include "predNA/EdgeMoveClusters.h"
#include "tools/InputParser.h"
#include "tools/StringTool.h"
#include "predNA/NuEnergyCalculator.h"
#include "predNA/EdgeClusterRegister.h"

namespace NSPpredNA {

class NuNode;
class EdgePartition;
class NuEdge;
class NuTree;
class SamplingGraph;
class NuGraph;
class graphInfo;


using namespace std;
using namespace NSPmodel;

class NuNode {

public:

	NuGraph* graph;
	int seqID;
	int baseType;

	bool connectToNeighbor;

	BaseConformer* baseConf;
	BaseConformer* baseConfTmp;

	RiboseConformer* riboseConf;
	RiboseConformer* riboseConfTmp;

	PhosphateConformer* phoConf;
	PhosphateConformer* phoConfTmp;

	double ene;
	double eneTmp;

	BaseConformerCG* baseConfCG;
	BaseConformerCG* baseConfCGTmp;

	RiboseConformerCG* riboseConfCG;
	RiboseConformerCG* riboseConfCGTmp;

	double bbcg;
	double bbcgTmp;

	double eneCG;
	double eneCGTmp;

	vector<int> neighborList;
	vector<int> connectionBreakPoints;

	vector<NuNode*> baseGroupA; //base coordinate not changed
	vector<NuNode*> riboGroupA; //ribose coordinate not changed

	vector<NuNode*> phoGroupA; //pho coordinate not changed
	vector<NuNode*> phoGroupC; //pho coordinate rotamer changed

	double samplingFreq;

	NuNode(int id, int baseType,LocalFrame& cs1, BaseRotamer* baseRot, RiboseRotamer* riboRot, AtomLib* atLib);
	NuNode(int id, int baseType,LocalFrame& cs1, BaseRotamerCG* baseRot, RiboseRotamerCG* riboRot, AtomLib* atLib);
	NuNode(int id, int baseType,LocalFrame& cs1, BaseRotamer* baseRot, BaseRotamerCG* baseRotCG, RiboseRotamer* riboRot, RiboseRotamerCG* riboRotCG, AtomLib* atLib);

	void updateEnergy(double clashRescale, double connectRescale);
	void updateEnergyCG(double clashRescale, double connectRescale);
	void updateNodeInformation(NuGraph* graph);
	void printNodeInfo();

	void updateRiboseRotamer(RiboseRotamer* rot, double clashRescale, double connectRescale);
	void acceptRotMutation();
	void clearRotMutation();
	void updateCoordinate(LocalFrame& cs);
	void acceptCoordMove();
	void clearCoordMove();
	double rotMutEnergy(double connectRescale);
	bool checkEnergy(double clashRescale, double connectRescale);

	void updateRiboseRotamerCG(RiboseRotamerCG* rot, double clashRescale, double connectRescale);
	void acceptRotMutationCG();
	void clearRotMutationCG();
	void updateCoordinateCG(LocalFrame& cs);
	void acceptCoordMoveCG();
	void clearCoordMoveCG();
	double rotMutEnergyCG();
	bool checkEnergyCG(double clashRescale, double connectRescale);
	int hbondNumTo(NuNode* other, AtomLib* atLib);

	vector<Atom*> toAtomList(AtomLib* atLib);
	vector<Atom*> toPhoAtomList(AtomLib* atLib);
	vector<Atom*> toAtomListOnlyBase(AtomLib* atLib);
	vector<Atom*> toAtomListWithPho(AtomLib* atLib);
	vector<Atom*> toAtomListCG(AtomLib* atLib);

	virtual ~NuNode();
};


class NuEdge {
public:

	NuGraph* graph;
	//the NuEdge is directed

	int indexA;
	int indexB;

	NuNode* nodeA;
	NuNode* nodeB;

	BasePairLib* pairLib;
	EdgeClusterRegister* ecr; //record the clusterIDs that occurred during sampling
	vector<EdgePartition*> epList; //random subTree starting from the edge
	EdgePartition* currentPartition;

	CsMove cm;
	CsMove cmTmp;
	int sep;

	double weight; //weight for MST
	double weightRand;

	EdgeMoveClusters* moveClusters;
	MixedNuPairCluster* moveSet;

	/*
	 * BB, BR, BP, RB, RR, RP, PB, PR, PP
	 *  0   1   2   3   4   5   6   7   8
	 */
 	double pairEne[9];
	double pairEneTmp[9]; //BB, BR, BP, RB, RR, RP, PB, PR, PP

	double pairEneCG[4]; //BB, BR, RB, RR
	double pairEneCGTmp[4]; //BB, BR, RB, RR

	double samplingFreq;

	NuEdge(NuNode* nodeA, NuNode* nodeB, NuGraph* graph);
	NuEdge(NuNode* nodeA, NuNode* nodeB, int sep, BasePairLib* pairLib, NuPairMoveSetLibrary* moveLib);

	void initNearNativeMoveClusters(double distanceCutoff=1.2);
	void initMoveClusters(BaseDistanceMatrix& dm, double distanceCutoff);
	void fixTemplateMove();

	void updateMoveClusters(double pContact=0.2); //according EdgeClusterRegister, add strong contact into move cluster

	void updateEnergy(double clashRescale, double connectRescale);
	void updateEnergyCG(double clashRescale, double connectRescale);

	void updateCsMove(CsMove& cm, double clashRescale, double connectRescale);

	void setRandomPartition();
	void recordLowEnergyCluster();

	double mutEnergy();
	
	void acceptMutation();
	void clearMutation();
	bool checkEnergy(double clashRescale, double connectRescale);
	bool checkReversePair();
	bool isWC();

	void updateCsMoveCG(CsMove& cm, double clashRescale, double connectRescale);
	double mutEnergyCG();
	void printMutEnergyCG();
	void acceptMutationCG();
	void clearMutationCG();
	bool checkEnergyCG(double clashRescale, double connectRescale);
	bool checkReversePairCG();
	void printPartition();
	virtual ~NuEdge();
};

class EdgePartition{

	public:
	/*
	 * edge indexA-indexB divide the NuTree into two trees, treeA and treeB
	 * we use BFS algorithm to traverses the tree and store all nodes and edges
	 * nodes of treeB are movable
	 */

	NuGraph* graph;
	NuEdge* targetEdge;

	vector<NuNode*> nodeListA; //coordinate not changed
	vector<NuEdge*> edgeListA;

	vector<NuNode*> nodeListB;
	vector<NuEdge*> edgeListB;

	vector<NuEdge*> edgeListX; //connection between nodeListA and nodeListB in SamplingGraph

	//positions where phosphate group need to be build
	vector<int> connectionBreakPoints;

	vector<NuNode*> phoGroupA; //pho coordinate not changed
	vector<NuNode*> phoGroupB; //pho coordinate changed, but rotamer not changed
	vector<NuNode*> phoGroupC; //pho coordinate rotamer changed

	string key;

	EdgePartition(NuEdge* edge, NuGraph* graph){
		this->targetEdge = edge;
		this->graph = graph;
	}

	void updatePartition(NuTree* tree);
	void updateEdgeListX(SamplingGraph* sg);

	void printPartition(){
		cout << "listA: ";
		for(int i=0;i<nodeListA.size();i++){
			cout << " " << nodeListA[i]->seqID ;
		}
		cout << endl;
		cout << "listB: ";
		for(int i=0;i<nodeListB.size();i++){
			cout << " " << nodeListB[i]->seqID ;
		}
		cout << endl;

		cout << "edgeListB: " << endl;
		for(int i=0;i<edgeListB.size();i++){
			cout << edgeListB[i]->nodeA->seqID << " - " << edgeListB[i]->nodeB->seqID << endl;
		}

		cout << "edgeListX: " << endl;
		for(int i=0;i<edgeListX.size();i++){
			cout << edgeListX[i]->nodeA->seqID << " - " << edgeListX[i]->nodeB->seqID << endl;
		}	
	}
};

class NuTree {
public:

	NuGraph* graph;

	bool* adjMtx; //adjacency matrix, N*N matrix, (N-1)*2 true points
	vector<NuEdge*> geList; 

	NuTree(NuGraph* graph);

//	void updateEdgeInfo(double clashRescale, double connectRescale);
//	void updateEdgeInfoCG(double clashRescale, double connectRescale);

	void printEdges();
	void printEdgeInfo(const string& output);
	void printEdgeInfo();
	void printTreeInfomation(const string& output); //print PDB in tre format

	virtual ~NuTree();
};

class SamplingGraph {

public:

	NuGraph* graph;

	vector<NuEdge*> geList;
	bool* adjMtx;
	vector<EdgePartition*> epList;


	int poolSize = 100000;
	double sampFreqNode;
	int randPoolNode[100000];
	double sampFreqEdge;
	int randPoolEdge[100000];
	double totalSamp;

	SamplingGraph(NuGraph* graph);


	void updatePartitionInfo();
	void updateSamplingInfo();
	void randomInit(double clashRescale, double connectRescale);
	void randomInitCG(double clashRescale, double connectRescale);

	void checkCsMove(){
		for(int i=0;i<geList.size();i++){
			NuEdge* edge = geList[i];
			LocalFrame csA = edge->nodeA->baseConf->cs1;
			LocalFrame csB = edge->nodeB->baseConf->cs1;
			LocalFrame csATmp = edge->nodeA->baseConfTmp->cs1;
			LocalFrame csBTmp = edge->nodeB->baseConfTmp->cs1;

			LocalFrame csB1 = csA + edge->cm;
			LocalFrame csBTmp1 = csATmp + edge->cmTmp;

			if(csB.origin_.distance(csB1.origin_) > 0.0001) {
				cout << "ERROR: edge " << edge->indexA << " " << edge->indexB << " move error" << endl; 
			}

			if(csBTmp.origin_.distance(csBTmp1.origin_) > 0.0001) {
				cout << "ERROR: edge " << edge->indexA << " " << edge->indexB << " tmp move error" << endl; 
			}
		}
	}

	void printNodeInfo();
	void printEdges();
	void printEdgeInfo(const string& output);
	void printEdgeInfo();
	void printTreeInfomation(const string& output); //print PDB in tre format

	graphInfo* runAtomicMC();
	void runCoarseGrainedMC(const string& output);
	virtual ~SamplingGraph();
};

class graphInfo{
public:
	int seqLen;
	int* seq;
	bool* connectToDownstream;
	int* sepTable;
	bool* fixed;
	NuNode** nodes;
	double ene;
	double nbEne;
	double nnbEne;
	double rms;
	AtomLib* atLib;

	//graphInfo(int seqLen, int* seq, bool* con, bool* fixed, NuNode** nodes, double ene, AtomLib* atLib);
	graphInfo(int seqLen, int* seq, bool* con, bool* fixed, NuNode** nodes, double ene, AtomLib* atLib, int mode); 

	void setRMS(double rms){
		this->rms = rms;
	}
	double rmsd(graphInfo* other);
	double rmsd(graphInfo* other, int pos);
	double rmsdCG(graphInfo* other);
	void printPDB(const string& outputFile);
	void printPDBWithPairMtx(const string& outputFile, BasePairLib* bpLib);
	void printAlignedPDB(graphInfo* alignTarget, const string& outputFile);
	void printDetailEnergy(const string& outputFile, BasePairLib* bpLib, AtomLib* atLib, RnaEnergyTable* et);
	void setNbEnergy(double e){
		this->nbEne = e;
	}
	void setNnbEnergy(double e){
		this->nnbEne = e;
	}
	virtual ~graphInfo();

};

class NuGraph {

public:
	int seqLen; //L
	int* seq; //sequenceType: 0~7
	int* wcPairPosID;
	int* stemIndex;

	bool* connectToDownstream; //bonded to 5' residue
	int* sepTable; //sequence seperation: -1, 0, 1, 2
	bool* fixed; //fixed residues, do not used for rms calculation

	vector<BaseRotamer*> initBaseRotList;
	vector<RiboseRotamer*> initRiboseRotList;
	vector<BaseRotamerCG*> initBaseRotCGList;
	vector<RiboseRotamerCG*> initRiboseRotCGList;

	NuNode** allNodes; //L nodes
	NuEdge** allEdges; //L*L edges
	vector<NuEdge*> geList; //L*(L-1)/2 edges
	RotamerLib* rotLib;
	AtomLib* atLib;

	BasePairLib* pairLib;
	EdgeMoveClustersLib* emcLib;
	NuPairMoveSetLibrary* moveLib;
	OrientationIndex* oi;
	RnaEnergyTable* et;
	graphInfo* initInfo;

	NuGraph(const string& inputFile, RotamerLib* rotLib, AtomLib* atLib, BasePairLib* pairLib, NuPairMoveSetLibrary* moveLib, EdgeMoveClustersLib* emcLib,  RnaEnergyTable* et);
	NuGraph(const string& inputFile, RotamerLib* rotLib, AtomLib* atLib, BasePairLib* pairLib);
	NuGraph(const string& inputFile, RotamerLib* rotLib, AtomLib* atLib, BasePairLib* pairLib, RnaEnergyTable* et, int InitMode);
	NuGraph(RNAPDB* pdb, RotamerLib* rotLib, AtomLib* atLib, BasePairLib* pairLib, NuPairMoveSetLibrary* moveLib, EdgeMoveClustersLib* emcLib,  RnaEnergyTable* et, const string& cnt);
	
	void init(const string& task, RNAPDB* pdb, const string& baseSeq, const string& baseSec, const string& csn, const string& cst, const string& cnt, const string& contactKey, vector<string>& ctList);

	void initPho();
	void initPho(PO3Builder* pb);
	void initForMC(const string& inputFile);
	void initForCGMC(const string& inputFile);
	void initForMST(const string& inputFile);
	void initForSingleResiduePrediction(const string& inputFile, int pos);
	void initRandWeight();
	void resetEdgeMoveToCurrentCluster();
	void MST_kruskal(NuTree* output, NuEdge* fixedEdge=NULL);
	void generateRandomEdgePartition(int roundNum=2);



	void updateEdgeMoveClusters(); //update move in cluster level
	void updateEdgeMoveSet(); //update move in GRID level


	void printAllEdge();

	void updateEnergy(double clashRescale, double connectRescale);
	void updateEnergyCG(double clashRescale, double connectRescale);
	void nodeListToPDBWithoutPho(vector<NuNode*> nodeList, RNAPDB* outpdb);

	void generateSubGraph(const string& inputFile, int corePos, int* subGraphPosList, int* fixedPositions, NuGraph* subGraph, vector<int>& outsubGraphPosList, vector<int>& outupdatePosList, vector<vector<int>>& Bclusters, vector<vector<int>>& Cclusters);
	void clusterContacts(int* array, int array_size, vector<vector<int>>& clusters);
	void dfsCluster(int* array, int array_size, int currentIdx, vector<bool>& visited, vector<int>& cluster) ;
	void mergeClusters(vector<vector<int>>& clusters, double threshold);
	double findMinDistanceBetweenClusters(const vector<int>& clusterA, const vector<int>& clusterB);
	string indexToBractString(int* index, int indexCount);
	void extractPDBAtoms(const string& inputFile, const string& outputFile, int* subGraphPosList, int subGraphCount);
	double minSquareDistance(int Node1,int Node2);

	string toContactMapHashKeyCG();

	void keyToContactMatrix(const string& key);
	void keyToMatrixFile(const string& key, const string& outfile);
	void printContactMatrix(const string& outfile);
	double keyAccuracy(const string& key);

	void checkCoordinate(){
		for(int i=0;i<seqLen;i++){
			if(this->allNodes[i]->baseConf->cs1.origin_.distance(this->allNodes[i]->baseConfTmp->cs1.origin_) > 0.01){
				cout << "node: " << i << " coord changed" << endl;
			}
		}
	}

	void checkEnergy(double clashRescale, double connectRescale);
	void checkEnergyCG(double clashRescale, double connectRescale);
	double totalEnergy(double clashRescale, double connectRescale);
	double nbEnergy(double clashRescale, double connectRescale);

	double nnbEnergy(double clashRescale, double connectRescale);
	double totalEnergyCG(double clashRescale, double connectRescale);
	double totalEnergyCGTmp(double clashRescale, double connectRescale);
	double totalEnergyTmp(double clashRescale, double connectRescale);
	double totalEnergy2();

	void printEnergy();
	void printEnergyCG(double clashRescale);
	void printBaseEnergyList(const string& outfile);
	void printPairwiseEnergy(const string& outfile);
	void printEdgeClusterRegister();
	void cgToAllAtom();

	double getTotalEnergyForModelSelection();

	graphInfo* getGraphInfo();
	graphInfo* getGraphInfoCG();
	graphInfo* getGraphInfo(double ene);
	graphInfo* getGraphInfoCG(double ene);

	virtual ~NuGraph();
};

} /* namespace NSPpredNA */

#endif /* PREDNA_NUGRAPH_H_ */
