/*
 * TestMCSimPara.cpp
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
using namespace NSPpredNA;
using namespace std;

int main(int argc, char** argv){
	clock_t start = clock();
	string inputFile = string(argv[1]);
	string paraFile = string(argv[2]);
	string output = string(argv[3]);

	RnaEnergyTable* et = new RnaEnergyTable(paraFile);

	BRFoldingTree* ft = new BRFoldingTree(inputFile, et);

	MCRun mc(ft);
}

