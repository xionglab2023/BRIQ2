/*
 * rnaGraph.cpp
 *
 *  Created on: 2020Äê11ÔÂ24ÈÕ
 *      Author: pengx
 */


#include "motif/RNAGraph.h"
#include <stdlib.h>
#include <iostream>
#include <time.h>

using namespace NSPmotif;
using namespace std;

int main(int argc, char** argv){

	string pdbFile = string(argv[1]);

	srand(time(0));

	RNAGraph* rg = new RNAGraph(pdbFile);
	rg->updateInitialEdge();
	rg->findHelix();
	rg->meltHelix();
	rg->meltLoop();
	rg->generateInitialClusters();


	vector<SubRNAGraph*> sgList = rg->findAllGraphs();

	for(int i=0;i<sgList.size();i++){
		sgList[i]->updateInfo();
		printf("B=%-2d F=%d Score=%7.3f\n", sgList[i]->B, sgList[i]->F, rg->getGraphScore(sgList[i]));
		sgList[i]->printFragmentString();
	}
	delete rg;
}

