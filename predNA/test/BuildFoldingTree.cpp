/*
 * BuildFoldingTree.cpp
 *
 *  Created on: 2022Äê5ÔÂ17ÈÕ
 *      Author: pengx
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

void printHelp(){
	cout << "Usage: buildFoldingTree $INPUT_FILE" << endl;
}

int main(int argc, char** argv){

	//Usage: briq_Predict $INPUT_FILE $OUTPUT_PDB
	if(argc != 2 || argv[1] == "-h")
	{
		printHelp();
		exit(0);
	}

	string inputFile = string(argv[1]);

	cout << "init folding tree" << endl;
	BRFoldingTree* ft = new BRFoldingTree(inputFile);

	cout << "print connection" << endl;
	ft->printConnections();

	delete ft;

}

