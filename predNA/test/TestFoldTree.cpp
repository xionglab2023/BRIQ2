/*
 * TestFoldTree.cpp
 *
 */


#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>

#include "model/StructureModel.h"
#include "predNA/BRFoldingTree.h"
#include "predNA/MCRun.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredna;
using namespace std;


int main(int argc, char** argv){


	clock_t start = clock();
	//Usage: briq_Predict $INPUT_FILE $OUTPUT_PDB

	string inputFile = string(argv[1]);
	string outputPDB = string(argv[2]);
	int randseed = atoi(argv[3]);

	srand(randseed);

	cout << "init energy table:" << endl;
	RnaEnergyTable* et = new RnaEnergyTable();

	clock_t end1 = clock();
	cout << "time1: " << (float)(end1-start)/CLOCKS_PER_SEC << "s" << endl;


	cout << "init folding tree" << endl;
	BRFoldingTree* ft = new BRFoldingTree(inputFile, et);

	clock_t end2 = clock();
	cout << "time2: " << (float)(end2-start)/CLOCKS_PER_SEC << "s" << endl;

	cout << "print connection" << endl;
	ft->printConnections();

	cout << "get tree info" << endl;
	BRTreeInfo* info = ft->getTreeInfo();
	cout << "run mc: " << endl;

	MCRun mc(ft);
	mc.simpleMC(outputPDB, false);

	delete et;
	delete ft;


}
