/*
 * GenerateDecoysRandInit.cpp
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
	RnaEnergyTable* et = new RnaEnergyTable(string(argv[2]));
	string output = string(argv[3]);
	int randSeed = atoi(argv[4]);


	srand(randSeed);
	BRFoldingTree* ft = new BRFoldingTree(inputFile, et);
	cout << "init mc run: "<< endl;
	MCRun mc(ft);



	mc.generateDecoysRandInit(output);

}



