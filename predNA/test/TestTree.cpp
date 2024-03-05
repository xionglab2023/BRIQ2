/*
 * TestTree.cpp
 *
 *  Created on: 2024,1,21
 *      Author: pengx
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


	BasePairLib* pairLib = new BasePairLib();
	RotamerLib* rotLib = new RotamerLib();
	AtomLib* atLib = new AtomLib();

	cout << "init graph" << endl;

    NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib);
	graph->initRandWeight();
	NuTree* tree = new NuTree(graph);
	graph->MST_kruskal(tree);

	cout << "tree: " << endl;
	tree->printEdgeInfo(output);


	delete pairLib;
	delete rotLib;
	delete atLib;
	delete tree;
	delete graph;


	clock_t end1 = clock();
	cout << "time1: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;

}


