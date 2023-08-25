/*
 * TestAlignment.cpp
 *
 *  Created on: 2020Äê12ÔÂ3ÈÕ
 *      Author: pengx
 */


#include "motif/MotifAlign.h"
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

	MotifAlign align(pdb1, pdb2);
	align.generateAlign();
	align.printResult();


}


