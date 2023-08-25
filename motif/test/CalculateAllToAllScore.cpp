/*
 * CalculateAllToAllScore.cpp
 *
 *  Created on: 2020Äê12ÔÂ10ÈÕ
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

	int indexA = atoi(argv[1]);

	ifstream file;

	string list = "/lustre/home/pengx/motif/mpdb/list";
	file.open(list.c_str(), ios::in);
	vector<RNAPDB*> pdbList;

	string s;
	while(getline(file, s)){
		RNAPDB* pdb = new RNAPDB("/lustre/home/pengx/motif/mpdb/" + s, "XXXX");
		pdbList.push_back(pdb);
	}
	file.close();

	int n = pdbList.size();

	ofstream out;
	char xx[200];
	sprintf(xx, "/lustre/home/pengx/motif/rmScore/%d.score", indexA);
	out.open(xx, ios::out);

	for(int i=indexA+1;i<n;i++){
		MotifAlign align(pdbList[indexA], pdbList[i]);
		align.generateAlign();
		double rmScore = align.getRMScore();
		out << indexA << " " << i << " " << rmScore << endl;
	}
	out.close();

}


