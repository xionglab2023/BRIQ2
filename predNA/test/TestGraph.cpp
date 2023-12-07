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
	//Usage: briq_Predict $INPUT_FILE $OUTPUT_PDB

	string inputFile = string(argv[1]);
	string outputPDB = string(argv[2]);
	int randseed = atoi(argv[3]);

	srand(randseed);

	/*
	BasePairLib* pairLib = new BasePairLib();
	NuPairMoveSetLibrary* moveLib = NULL;

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
	*/

	cout << "init graph" << endl;
	NuGraph* graph = new NuGraph(inputFile);
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

	cout << "run MC" << endl;
	tree->runAtomicMC();

	clock_t end1 = clock();
	cout << "time1: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;

}


