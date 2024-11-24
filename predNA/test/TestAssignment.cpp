/*
 * TestAssignment.cpp
 *
 */

#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include "model/AssignRNASS.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>

#include "model/StructureModel.h"
#include "predNA/BRFoldingTree.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredNA;
using namespace std;

int main(int argc, char** argv){
	/*
	 * input: $PDBFILE
	 */

	if(argc != 3 || argv[1] == "-h")
	{
		cout << "Usage: BRiQ_assignSS $PDBFILE $OUTFILE" << endl;
		exit(0);
	}

	//cout << "start" << endl;
	string pdbFile = string(argv[1]);
	string outFile = string(argv[2]);

	RNAPDB* pdb = new RNAPDB(pdbFile, "xxxx");
	AtomLib* atLib = new AtomLib();

	//cout << "init ar" << endl;
	AssignRNASS* ar = new AssignRNASS(pdb, atLib);
	//cout << "assign " << endl;
	ar->printInfo(outFile);

	delete pdb;
	delete atLib;
	delete ar;

}


