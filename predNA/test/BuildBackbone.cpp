/*
 * BuidBackbone.cpp
 *
 */


#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "model/StructureModel.h"
#include "predNA/BackboneModelingTemplate.h"
#include "predNA/MCRun.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredna;
using namespace std;

int main(int argc, char** argv){

	string pdbFile = string(argv[1]);
	//string paraFile = string(argv[2]);
	string output = string(argv[2]);

	cout << "start: " << endl;


	cout << "init bm" << endl;
	BackboneModelingTemplate bm(pdbFile);

	cout << "runMC" << endl;
	//bm.printEnergyDetail();
	//bm.debug();
	bm.runMC();

	bm.printEnergyDetail();

	//bm.printDiheds(output);

	BRTreeInfo* info = bm.toTreeInfo();
	info->printPDB(output);
}

