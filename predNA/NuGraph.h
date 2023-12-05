/*
 * NuGraph.h
 *
 *  Created on: 2023Äê11ÔÂ15ÈÕ
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
#include "model/StructureModel.h"
#include "model/AtomLib.h"
#include "model/RotamerLib.h"
#include "model/BaseRotamer.h"
#include "model/RiboseRotamer.h"
#include "model/PhosphateRotamer.h"
#include "model/BasePairLib.h"
#include "predNA/NuMoveSet.h"
#include "predNA/EdgeInformation.h"
#include "tools/InputParser.h"
#include "tools/StringTool.h"


namespace NSPpredNA {

class NuNode;
class NuEdge;
class EdgeInformation;
class NuTree;
class NuGraph;


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

	double eneCG;
	double eneCGTmp;

	vector<int> neighborList;
	vector<int> connectionBreakPoints;

	NuNode(int id, int baseType,LocalFrame& cs1, RiboseRotamer* riboRot, AtomLib* atLib);

	void updateNodeInformation(NuTree* tree);
	void updateRiboseRotamer(RiboseRotamer* rot);
	void acceptRotMutation();
	void clearRotMutation();
	void updateCoordinate(LocalFrame& cs);
	void acceptCoordMove();
	void clearCoordMove();
	double mutEnergy();
	void updateRiboseRotamerCG(RiboseRotamerCG* rot);
	void acceptRotMutationCG();
	void clearRotMutationCG();
	void updateCoordinateCG(LocalFrame& cs);
	void acceptCoordMoveCG();
	void clearCoordMoveCG();
	double mutEnergyCG();


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

	double weight; //weight for MST
	double weightRand;


	EdgeInformation* ei;
	MixedNuPairCluster* moveSet;

	NuEdge(NuNode* nodeA, NuNode* nodeB, NuGraph* graph);

	NuEdge(NuNode* nodeA, NuNode* nodeB, int sep, BasePairLib* pairLib, NuPairMoveSetLibrary* moveLib);

	NuEdge& operator=(const NuEdge& other){
		this->graph = other.graph;
		this->indexA = other.indexA;
		this->indexB = other.indexB;
		this->nodeA = other.nodeA;
		this->nodeB = other.nodeB;
		this->cm = other.cm;
		this->cmTmp = other.cmTmp;
		this->sep = other.sep;

		this->weight = other.weight;
		this->weightRand = other.weightRand;
		this->ei = other.ei;
		this->moveSet = other.moveSet;

		for(int i=0;i<other.nodeListA.size();i++){
			this->nodeListA.push_back(other.nodeListA[i]);
		}
		for(int i=0;i<other.nodeListB.size();i++){
			this->nodeListB.push_back(other.nodeListB[i]);
		}
		for(int i=0;i<other.edgeListA.size();i++){
			this->edgeListA.push_back(other.edgeListA[i]);
		}
		for(int i=0;i<other.edgeListB.size();i++){
			this->edgeListB.push_back(other.edgeListB[i]);
		}
		for(int i=0;i<other.connectionBreakPoints.size();i++){
			this->connectionBreakPoints.push_back(other.connectionBreakPoints[i]);
		}
		return *this;
	}

	bool operator<(const NuEdge& other){
		return this->weightRand  < other.weightRand;
	}

	void initNativeMoveSet();
	void fixNaiveMove();

	void updateSubTrees(NuTree* tree);
	void updateCsMove(CsMove& cm);
	double mutEnergy();
	void acceptMutation();
	void clearMutation();

	void updateCsMoveCG(CsMove& cm);
	double mutEnergyCG();
	void acceptMutationCG();
	void clearMutationCG();

	virtual ~NuEdge();
};



class NuTree {
public:

	NuGraph* graph;
	bool* adjMtx; //adjacency matrix, L*L matrix, (L-1)*2 true points
	vector<NuEdge*> geList; //edge list, (L-1)
	NuTree(NuGraph* graph);
	void printEdges();
	virtual ~NuTree();
};

class NuGraph {

public:
	int len; //L
	int* seq; //sequenceType: 0~7
	int* wcPairPosID;
	bool* connectToDownstream;
	int* sepTable; //sequence seperation: -1, 0, 1, 2

	vector<RiboseRotamer*> initRiboseRotList;

	NuNode** allNodes; //L nodes
	NuEdge** allEdges; //L*L edges

	vector<NuEdge*> geList; //connected edges

	RotamerLib* rotLib;
	AtomLib* atLib;
	NuPairMoveSetLibrary* moveLib;
	BasePairLib* pairLib;

	NuGraph(const string& inputFile);

	void init(const string& inputFile);
	void initRandWeight();
	void MST_kruskal(NuTree* output);

	void printAllEdge();

	virtual ~NuGraph();
};

} /* namespace NSPpredNA */

#endif /* PREDNA_NUGRAPH_H_ */
