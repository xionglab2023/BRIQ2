/*
 * TestAlignScore.cpp
 *
 *  Created on: 2020Äê12ÔÂ9ÈÕ
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
	srand(atoi(argv[1]));
	//double d0 = atof(argv[2]);

	ifstream file;

	string list = "/lustre/home/pengx/motif/mpdb/list";
	file.open(list.c_str(), ios::in);
	vector<RNAPDB*> pdbList;

	string s;
	while(getline(file, s)){
		RNAPDB* pdb = new RNAPDB("/lustre/home/pengx/motif/mpdb/" + s, "XXXX");
		if(pdb->getBaseList().size() > 50) {
			cout << s << endl;
			delete pdb;
			continue;
		}
		pdbList.push_back(pdb);
	}
	file.close();


	int n = pdbList.size();
	for(int i=0;i<10000;i++){
		int idA = rand()%n;
		int idB = rand()%n;
		int lenSum = pdbList[idA]->getBaseList().size() + pdbList[idB]->getBaseList().size();


		for(double d0=0.02; d0<=1.4; d0+=0.02){
			if(d0 <= 0) continue;
			MotifAlign align(pdbList[idA], pdbList[idB]);
			align.generateAlign();
			align.printResult();
		}
	}



}

