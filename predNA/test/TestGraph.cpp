/*
 * TestGraph.cpp
 *
 *  Created on: 2023Äê12ÔÂ1ÈÕ
 *      Author: nuc
 */



#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>

#include "model/StructureModel.h"
#include "predNA/NuGraph.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredNA;
using namespace std;


int main(int argc, char** argv){


	clock_t start = clock();

	string inputFile = string(argv[1]);
	string output = string(argv[2]);
	int randseed = atoi(argv[3]);

	srand(randseed);

	/*
	BasePairLib* pairLib = new BasePairLib();
	NuPairMoveSetLibrary* moveLib = NULL;
	RnaEnergyTable* et = new RnaEnergyTable();

	AtomLib* atLib = new AtomLib();
	RNAPDB* pdb = new RNAPDB("/public/home/pengx/cpp/briqx/demo/gcaa/gcaa.pdb", "xxx");
	vector<RNABase*> baseList = pdb->getBaseList();
	vector<NuNode*> nodeList;
	vector<NuEdge*> edgeList;

	for(int i=0;i<baseList.size();i++)
	{
		cout << "init node" << i << endl;
		RNABase* base = baseList[i];
		LocalFrame cs = base->getCoordSystem();
		RiboseRotamer* riboRot = new RiboseRotamer(base);
		NuNode* node = new NuNode(i, base->baseTypeInt, cs, riboRot, atLib);
		nodeList.push_back(node);
	}

	for(int i=0;i<baseList.size();i++){
		for(int j=0;j<baseList.size();j++){

			int sep = j-i;
			if(abs(sep) > 1) sep = 2;
			cout << "init edge: " << i << " " << j << endl;
			NuEdge* edge = new NuEdge(nodeList[i], nodeList[j], sep, pairLib, moveLib);
			cout << "move set" << endl;
			edge->initNativeMoveSet();
			edgeList.push_back(edge);
		}
	}

	for(int i=0;i<baseList.size();i++){
		BaseConformer* confA = nodeList[i]->baseConf;
		for(int j=i+2;j<baseList.size();j++){
			BaseConformer* confB = nodeList[j]->baseConf;
			double e1 = nuBaseBaseEnergy(confA, confB, 2, et);
			double e2 = nuBaseBaseEnergy(confB, confA, 2, et);
			if(abs(e1+e2) > 0){
				printf("baseA: %d baseB: %d e1: %7.3f e2: %7.3f\n", i, j, e1, e2);
			}
		}
	}
	*/


	//NuPairMoveSetLibrary* moveSet = new NuPairMoveSetLibrary();


	/*
	int typeA = 0;
	int typeB = 2;
	int sep = 2;
	int pairType = typeA*4+typeB;
	int clusterID = 2;

	OrientationIndex* oi = new OrientationIndex();
	BasePairLib* pairLib = new BasePairLib();

	EdgeInformation* ei = new EdgeInformation(sep, typeA, typeB, pairLib);
	ei->setUniqueCluster(2, pairLib);
	MixedNuPairCluster* mc = new MixedNuPairCluster(sep, typeA*4+typeB, moveSet);
	mc->updateEdgeInformation(ei);

	for(int i=0;i<100;i++){
		CsMove cm = mc->getRandomMove();
		BaseDistanceMatrix dm(cm);
		int clusterID = pairLib->getPairType(dm, typeA, typeB, sep);
		double d = pairLib->nnbDMClusterCenters[typeA*4+typeB][clusterID].distanceTo(dm);
		printf("%-3d %5.3f\n", clusterID, d);
	}
	delete moveSet;
	delete pairLib;
	delete ei;
	delete mc;
	*/

	/*
	NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib);
	NuTree* tree = new NuTree(graph);
	graph->MST_kruskal(tree);
	tree->printEdgeInfo(output);
	*/


	/*
	 * check tree connection
	*/

	BasePairLib* pairLib = new BasePairLib();
	RotamerLib* rotLib = new RotamerLib();
	AtomLib* atLib = new AtomLib();
	NuPairMoveSetLibrary* moveLib = new NuPairMoveSetLibrary();
	RnaEnergyTable* et = new RnaEnergyTable();


	cout << "init graph" << endl;
	NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, et);
	graph->initRandWeight();
	cout << "all edge" << endl;
	graph->printAllEdge();
	NuTree* tree = new NuTree(graph);
	graph->MST_kruskal(tree);
	cout << "tree: " << endl;
	tree->printEdges();
	cout << "adde node info" << endl;
	tree->updateNodeInfo();
	for(int i=0;i<graph->seqLen;i++){
		graph->allNodes[i]->printNodeInfo();
	}

	cout << "edgeInfo: " << endl;
	tree->updateEdgeInfo();
	for(int i=0;i<tree->geList.size();i++){
		cout << "edge: " << tree->geList[i]->indexA << "-" << tree->geList[i]->indexB << endl;
		tree->geList[i]->printPartition();
	}

	cout << "update sampling info" << endl;
	tree->updateSamplingInfo();
	tree->printNodeInfo();

	cout << "run MC" << endl;
	tree->runAtomicMC(output);

	delete pairLib;
	delete rotLib;
	delete atLib;
	delete moveLib;
	delete et;
	delete tree;
	delete graph;


	clock_t end1 = clock();
	cout << "time1: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;

}


