/*
 * runRnaXp.cpp
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

	//EnergyTable* et = NULL;
	RnaEnergyTable* et = new RnaEnergyTable(string(argv[2]));

	string outputRMS = string(argv[3]);
	string outputPDB = string(argv[4]);

	int randSeed = atoi(argv[5]);
	srand(randSeed);

	cout << "init folding tree" << endl;
	BRFoldingTree* ft = new BRFoldingTree(inputFile, et);
	BRTreeInfo* info = ft->getTreeInfo();

	cout << "run mc: " << endl;

	MCRun mc(ft);
	mc.simpleMC(outputPDB, false);

	info = ft->getTreeInfo();

	cout << "get rms: " << endl;
	ofstream out;
	out.open(outputRMS, ios::out);
	out << info->rmsd(mc.init) << endl;
	out << info->ene << endl;
	ft->printDetailEnergy(out);
	clock_t end = clock();
	out << "time: " << (float)(end-start)/CLOCKS_PER_SEC << "s" << endl;
	out.close();

	delete et;
	delete ft;

}


