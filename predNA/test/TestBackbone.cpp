/*
 * TestBackbone.cpp
 *
 */

#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>

#include "model/StructureModel.h"
#include "predNA/MCRun.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredna;
using namespace std;


int main(int argc, char** argv){


	clock_t start = clock();

	string inputFile = string(argv[1]);
	string pdbID = "XXXX";
	RNAPDB pdb(inputFile, pdbID);
	vector<RNABase*> baseList = pdb.getBaseList();

	RNABase* baseA = baseList[1];
	RNABase* baseB = baseList[2];

	RotamerLib* rotLib = new RotamerLib();
	RiboseRotamer* rotA = new RiboseRotamer(baseA);
	RiboseRotamer* rotB = new RiboseRotamer(baseB);

	PhosphateRotamer* phoRotA = new PhosphateRotamer(0.0, 0.0);
	PhosphateRotamer* phoRotB = new PhosphateRotamer(0.0, 0.0);


	BRNode* nodeA = new BRNode(baseA, rotA, phoRotA, rotLib);
	BRNode* nodeB = new BRNode(baseB, rotB, phoRotB, rotLib);

	LocalFrame csA1 = nodeA->baseConf->cs1;
	LocalFrame csA2 = nodeA->riboseConf->cs2;
	LocalFrame csA3 = nodeA->riboseConf->cs1;

	LocalFrame csB1 = nodeB->baseConf->cs1;
	LocalFrame csB2 = nodeB->riboseConf->cs2;
	LocalFrame csB3 = nodeB->riboseConf->cs1;

	CsMove cm = csB3 - csA2;

	CsMove mv = rotA->mv12 + cm + rotB->mv31;
	LocalFrame csBB = csA1 + mv;

	csB1.print();
	cout << endl;
	csBB.print();

}


