/*
 * optModel.cpp
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
	string outFile = string(argv[3]);
	string outPDB = string(argv[4]);

	srand(atoi(argv[5]));

	double t0 = 0.5;
	double kStep = 1.0;

	BRFoldingTree* ft = new BRFoldingTree(inputFile, et);
	MCRun mc(ft);
	mc.optimize(t0, kStep);

	clock_t end = clock();

	ofstream of;
	of.open(outFile.c_str(), ios::out);
	ft->printDetailEnergy(of);
	of << "time: " << (float)(end-start)/CLOCKS_PER_SEC << "s" << endl;
	of.close();


	BRTreeInfo* info = ft->getTreeInfo();
	info->printPDB(outPDB);



	//treeInfo* out = run.optFromRandomInit(logFile);
}

