/*
 * ConnectInfoToRMS.cpp
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
	cout << "init energy table:" << endl;
	string inputFile = string(argv[1]);
	RnaEnergyTable* et = new RnaEnergyTable();
	string keyfile = string(argv[2]);
	string outfile = string(argv[3]);

	BRFoldingTree* ft = new BRFoldingTree(inputFile, et);
	ft->keyInfoToRMS(keyfile, outfile);


}


