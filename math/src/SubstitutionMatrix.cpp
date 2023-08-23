/*
 * SubstitutionMatrix.cpp
 *
 *  Created on: 2021Äê12ÔÂ10ÈÕ
 *      Author: pengx
 */

#include "math/SubstitutionMatrix.h"

namespace NSPmath {

SubstitutionMatrix::SubstitutionMatrix() {
	// TODO Auto-generated constructor stub
	string path = NSPdataio::datapath();
	string subFile = path + "subMatrix/sub62.mtx";
	ifstream f;
	f.open(subFile.c_str(), ios::in);
	string line;
	vector<string> spt;
	int id = 0;
	while(getline(f, line)){
		NSPtools::splitString(line, " ", &spt);
		for(int j=0;j<20;j++){
			this->subM[id][j] = atof(spt[j].c_str());
		}
		id++;
	}
	f.close();
}

SubstitutionMatrix::~SubstitutionMatrix() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPmath */
