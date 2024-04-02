/*
 * NuGraph.h
 *
 *  Created on: 2023��11��15��
 *      Author: nuc
 */

#ifndef PREDNA_NUGRAPH_H_
#define PREDNA_NUGRAPH_H_

#include <vector>
#include <string>
#include <map>
#include <queue>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <memory>
#include "model/StructureModel.h"
#include "model/AtomLib.h"
#include "model/RotamerLib.h"
#include "model/BaseRotamer.h"
#include "model/RiboseRotamer.h"
#include "model/PhosphateRotamer.h"
#include "model/BasePairLib.h"
#include "geometry/RMSD.h"
#include "predNA/NuMoveSet.h"
#include "predNA/EdgeInformation.h"
#include "tools/InputParser.h"
#include "tools/StringTool.h"
#include "NuEnergyCalculator.h"

namespace NSPpredNA {

class NuNode;
class NuEdge;
class EdgeInformation;
class NuTree;
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

	vector<NuNode*> baseGroupA; //base coordinate not changed, without masked nodes
	vector<NuNode*> riboGroupA; //ribose coordinate not changed, without masked nodes

	vector<NuNode*> phoGroupA; //pho coordinate not changed, without masked nodes
	vector<NuNode*> phoGroupC; //pho coordinate rotamer changed

	double samplingFreq;

	NuNode(int id, int baseType,LocalFrame& cs1, BaseRotamer* baseRot, RiboseRotamer* riboRot, AtomLib* atLib);
	NuNode(int id, int baseType,LocalFrame& cs1, BaseRotamerCG* baseRot, RiboseRotamerCG* riboRot, AtomLib* atLib);
	NuNode(int id, int baseType,LocalFrame& cs1, BaseRotamer* baseRot, BaseRotamerCG* baseRotCG, RiboseRotamer* riboRot, RiboseRotamerCG* riboRotCG, AtomLib* atLib);

	void updateNodeInformation(NuTree* tree);
	void updateNodeInformationCG(NuTree* tree);
	void printNodeInfo();

	void updateRiboseRotamer(RiboseRotamer* rot);
	void acceptRotMutation();
	void clearRotMutation();
	void updateCoordinate(LocalFrame& cs);
	void acceptCoordMove();
	void clearCoordMove();
	double rotMutEnergy();
	bool checkEnergy();

	void updateRiboseRotamerCG(RiboseRotamerCG* rot);
	void acceptRotMutationCG();
	void clearRotMutationCG();
	void updateCoordinateCG(LocalFrame& cs);
	void acceptCoordMoveCG();
	void clearCoordMoveCG();
	double rotMutEnergyCG();
	bool checkEnergyCG();


	vector<Atom*> toAtomList(AtomLib* atLib);
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

	CsMove cm;
	CsMove cmTmp;
	int sep;

	/*
	 * edge indexA-indexB divide the NuTree into two trees, treeA and treeB
	 * we use BFS algorithm to traverses the tree and store all nodes and edges
	 * nodes of treeB are movable
	 */

	vector<NuNode*> nodeListA;
	vector<NuEdge*> edgeListA;

	vector<NuNode*> nodeListB;
	vector<NuEdge*> edgeListB;

	//positions where phosphate group need to be build
	vector<int> connectionBreakPoints;

	vector<NuNode*> phoGroupA; //pho coordinate not changed
	vector<NuNode*> phoGroupB; //pho coordinate changed, but rotamer not changed
	vector<NuNode*> phoGroupC; //pho coordinate rotamer changed

	double weight; //weight for MST
	double weightRand;

	EdgeInformation* ei;
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

	void initNearNativeMoveSet();
	void fixNaiveMove();

	void updateEdgeInfo(NuTree* tree);
	void updateEdgeInfoCG(NuTree* tree);

	void updateCsMove(CsMove& cm);
	double mutEnergy();
	void acceptMutation();
	void clearMutation();
	bool checkEnergy();
	bool checkReversePair();

	void updateCsMoveCG(CsMove& cm);
	double mutEnergyCG();
	void acceptMutationCG();
	void clearMutationCG();
	bool checkEnergyCG();
	bool checkReversePairCG();
	void printPartition();
	virtual ~NuEdge();
};

class NuTree {
public:

	NuGraph* graph;

	bool* masked; //
	bool* adjMtx; //adjacency matrix, N*N matrix, (N-1)*2 true points
	vector<NuEdge*> geList; //edge list, (N-1)

	int poolSize = 100000;
	double sampFreqNode;
	int randPoolNode[100000];
	double sampFreqEdge;
	int randPoolEdge[100000];
	double totalSamp;


	NuTree(NuGraph* graph);
	void updateNodeInfo();
	void updateNodeInfoCG();
	void printNodeInfo();
	void updateEdgeInfo();
	void updateEdgeInfoCG();
	void updateSamplingInfo();
	
	void randomInit();
	void randomInitCG();

	void printEdges();
	void printEdgeInfo(const string& output);
	void printEdgeInfo();

	graphInfo* runAtomicMC();
	void runCoarseGrainedMC(const string& output);
	virtual ~NuTree();
};

class graphInfo{
public:
	int seqLen;
	int* seq;
	bool* connectToDownstream;
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
	void printAlignedPDB(graphInfo* alignTarget, const string& outputFile);
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
	bool* connectToDownstream; //bonded to 5' residue
	int* sepTable; //sequence seperation: -1, 0, 1, 2
	bool* masked; //masked residues do not participate in energy calculation and structure sampling, do not exist in NuTree
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
	NuPairMoveSetLibrary* moveLib;
	RnaEnergyTable* et;
	graphInfo* initInfo;

	NuGraph(const string& inputFile, RotamerLib* rotLib, AtomLib* atLib, BasePairLib* pairLib, NuPairMoveSetLibrary* moveLib, RnaEnergyTable* et);
	NuGraph(const string& inputFile, RotamerLib* rotLib, AtomLib* atLib, BasePairLib* pairLib);

	void init(const string& task, const string& pdbFile, const string& baseSeq, const string& baseSec, const string& cst, const string& chainBreak);
	void initPho();
	void initPho(PO3Builder* pb);
	void initForMC(const string& inputFile);
	void initForCGMC(const string& inputFile);
	void initForMST(const string& inputFile);
	void initForSingleResiduePrediction(const string& inputFile, int pos);
	void initRandWeight();
	void MST_kruskal(NuTree* output);
	void printAllEdge();
	void checkEnergy();
	void checkEnergyCG();

	double totalEnergy();

	double nbEnergy();
	double nnbEnergy();

	double totalEnergyCG();
	double totalEnergyCGTmp();
	
	double totalEnergyTmp();
	double totalEnergy2();

	void printEnergy();
	void printEnergyCG();

	void cgToAllAtom();

	graphInfo* getGraphInfo();
	graphInfo* getGraphInfoCG();

	virtual ~NuGraph();
};

} /* namespace NSPpredNA */

#endif /* PREDNA_NUGRAPH_H_ */
