/*
 * RNAScoring.cpp
 *
 */


#include "geometry/localframe.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

#include "forcefield/RnaEnergyTable.h"
#include "model/StructureModel.h"
#include "predNA/BRFoldingTree.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredNA;
using namespace std;

int main(int argc, char** argv){
	string inputFile = string(argv[1]);
	cout << "init energy table:" << endl;
	ForceFieldPara* para = new ForceFieldPara();

	RnaEnergyTable* et = new RnaEnergyTable(para);

	cout << "init folding tree" << endl;
	BRFoldingTree* ft = new BRFoldingTree(inputFile, et);

	ft->printDetailEnergy();

	delete et;
	delete ft;

}



