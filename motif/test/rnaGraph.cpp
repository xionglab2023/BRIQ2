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
	string outFile = string(argv[2]);

	ofstream out;
	out.open(outFile.c_str(), ios::out);


	srand(time(0));

	set<string> results;

	for(int k=0;k<10;k++) {
		RNAGraph* rg = new RNAGraph(pdbFile);
		rg->updateInitialEdge();
		rg->findHelix();
		rg->meltHelix();
		rg->meltLoop();
		rg->generateInitialClusters();
		vector<SubRNAGraph*> sgList = rg->findAllGraphs();
		for(int i=0;i<sgList.size();i++){
			results.insert(sgList[i]->toString());
		}
		delete rg;
	}


	set<string>::iterator it;

	for(it = results.begin();it!= results.end();++it){
		out << *it << endl;
	}
	out.close();


}

