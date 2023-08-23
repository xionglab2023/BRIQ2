/*
 * ResScRotamerLibMini.cpp
 *
 *  Created on: 2022Äê7ÔÂ16ÈÕ
 *      Author: pengx
 */

#include "model/ResScRotamerLibMini.h"

namespace NSPmodel {

ResScRotamerLibMini::ResScRotamerLibMini() {
	// TODO Auto-generated constructor stub
	this->atLib = new AtomLib();
	string path = NSPdataio::datapath();
	ResName rn;
	ifstream file;
	string s;
	string fileName;
	vector<string> spt;
	int index;
	char xx[200];

	for(int i=0;i<20;i++){
		string tri = rn.intToTri(i);
		if(tri == "GLY") {
			this->rotList[5].push_back(new ResScRotamer());
			continue;
		}
		fileName = path + "scRotamer/rotMin/" + tri +".rot";
		file.open(fileName.c_str(), ios::in);
		if(! file.is_open()) {
			cout << "fail to open sidechain rotamer library file " << fileName << endl;
			exit(1);
		}

		getline(file, s);
		int index = 0;
		while(getline(file, s)){
			ResScRotamer* rot = new ResScRotamer(s, atLib);
			rot->rotID = index;
			this->rotList[i].push_back(rot);
			index++;
		}
		file.close();
	}

	for(int i=0;i<20;i++){
		this->rotNum.push_back(this->rotList[i].size());
	}

	int idA, idB;
	double ene;

	for(int i=0;i<20;i++){
		for(int j=0;j<1000;j++){
			for(int k=0;k<80;k++){
				eRot[i][j][k] = 0.0;
			}
		}
	}

	for(int aa=0;aa<20;aa++){
		string restype = rn.intToTri(aa);
		fileName = path + "scRotamer/eneMin/" + restype + ".ene";
		file.open(fileName.c_str(), ios::in);
		if(! file.is_open()) {
			cout << "fail to open sidechain rotamer energy file " << fileName << endl;
			exit(1);
		}
		while(file >> idA >> idB >> ene){
				this->eRot[aa][idA][idB] = ene;
		}
		file.close();
	}

}

ResScRotamerLibMini::~ResScRotamerLibMini() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPmodel */
