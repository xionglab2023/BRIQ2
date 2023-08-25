/*
 * TestSimpleAlign.cpp
 *
 *  Created on: 2023Äê8ÔÂ25ÈÕ
 *      Author: nuc
 */


#include "motif/SimpleAlign.h"
#include <stdlib.h>
#include <iostream>
#include <time.h>

using namespace NSPmodel;
using namespace NSPmotif;
using namespace std;

int main(int argc, char** argv){

	string pdbFile1 = string(argv[1]);
	string pdbFile2 = string(argv[2]);

	RNAPDB* pdb1 = new RNAPDB(pdbFile1, "xxxx");
	RNAPDB* pdb2 = new RNAPDB(pdbFile2, "xxxx");

	SimpleAlign align(pdb1, pdb2);
	double score = align.getRMScore();
	//printf("%5.3f\n", score);
	align.printResult();

}



