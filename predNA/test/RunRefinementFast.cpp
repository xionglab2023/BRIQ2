/*
 * RunRefinementFast.cpp
 *
 *  Created on: 2022Äê6ÔÂ6ÈÕ
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

void checkInputFile(const string& inputFile){
	ifstream file;
	cout << "check input file: " << endl;
	file.open(inputFile.c_str(), ios::in);
	if(!file.is_open()){
		cout << "fail to open inputFile: " << inputFile << endl;
		exit(0);
	}
	file.close();

	NSPtools::InputParser input(inputFile);
	string pdbFile = input.getValue("pdb");
	string baseSeq = input.getValue("seq");
	string baseSec = input.getValue("sec");
	string nwcSec = input.getValue("nwc");

	if(pdbFile == ""){
		cout << "invalid inputFile: initial pdb required" << endl;
		exit(0);
	}

	if(baseSeq == ""){
		cout << "invalid inputFile: base sequence required" << endl;
		exit(0);
	}

	if(baseSec == ""){
		cout << "invalid inputFile: Watson-Crick pairing sequence required" << endl;
		exit(0);
	}

	if(nwcSec == ""){
		cout << "invalid inputFile: non-Watson-Crick pairing sequence required" << endl;
		exit(0);
	}

	if(baseSec.length() != baseSeq.length()){
		cout << "the length of Watson-Crick pairing sequence must equal to the length of base sequence" << endl;
		exit(0);
	}

	if(nwcSec.length() != baseSeq.length()){
		cout << "the length of non-Watson-Crick pairing sequence must equal to the length of base sequence" << endl;
		exit(0);
	}


	file.open(pdbFile.c_str(), ios::in);
	if(!file.is_open()){
		cout << "fail to open pdb file: " << pdbFile << endl;
		exit(0);
	}
}

void printHelp(){
	cout << "Usage: briq_Refinement_fast $INPUT_FILE $OUTPUT_PDB $RANDOMSEED" << endl;
}

int main(int argc, char** argv){

	if(argc != 4 || argv[1] == "-h")
	{
		printHelp();
		exit(0);
	}

	string inputFile = string(argv[1]);
	string outPDB = string(argv[2]);
	int randseed = atoi(argv[3]);

	srand(randseed);

	checkInputFile(inputFile);

	cout << "init energy table:" << endl;
	RnaEnergyTable* et = new RnaEnergyTable();

	double t0 = 0.4;
	double kStep = 0.4;

	cout << "run refinement: " << endl;
	BRFoldingTree* ft = new BRFoldingTree(inputFile, et);
	MCRun mc(ft);
	mc.optimize(t0, kStep);

	BRTreeInfo* info = ft->getTreeInfo();
	info->printPDB(outPDB);

}





