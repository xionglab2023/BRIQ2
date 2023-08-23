/*
 * ResBBRotamerLib.cpp
 *
 */

#include "model/ResBBRotamerLib.h"

namespace NSPmodel {
ResBBRotamerLib::ResBBRotamerLib() {
	// TODO Auto-generated constructor stub
	string path = NSPdataio::datapath();
	ResName rn;

	this->atLib = new AtomLib();
	ifstream file;
	string s;

	string libFile = path + "/bbRotamer/rot1k/all.nearby";
	int id = 0;
	file.open(libFile.c_str(), ios::in);
	if(! file.is_open()){
		cout << "fail to open backbone rotamer library file " << libFile << endl;
		exit(1);
	}

	id = 0;
	vector<string> spt;
	while(getline(file, s)){
		splitString(s, " ", &spt);
		for(int j=0;j<50;j++){
			this->neighborLib1k[id][j] = atoi(spt[j].c_str());
		}
		id++;
	}
	file.close();


	libFile = path + "/bbRotamer/rot1k/all.Lv3Neighbor";
	file.open(libFile.c_str(), ios::in);
	if(! file.is_open()){
		cout << "fail to open backbone rotamer library file " << libFile << endl;
		exit(1);
	}
	id = 0;
	while(getline(file, s)){
		splitString(s, " ", &spt);
		for(int j=0;j<12;j++){
			this->neighborLib1w[id][j] = atoi(spt[j].c_str());
		}
		id++;
	}
	file.close();


	id = 0;
	libFile = path + "/bbRotamer/rot1k/all.rot";
	file.open(libFile.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open backbone rotamer library file " << libFile << endl;
		exit(1);
	}
	while(getline(file, s)){
		for(int aa=0;aa<20;aa++) {
			ResBBRotamer* rot = new ResBBRotamer(s, atLib);
			rot->setAAType(aa, atLib);
			this->allRotLib1k[aa][id] = rot;
		}
		id++;
	}
	file.close();


	libFile = path + "/bbRotamer/rot1w/all.rot";
	file.open(libFile.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open backbone rotamer library file " << libFile << endl;
		exit(1);
	}
	id = 0;
	while(getline(file, s)){
		for(int aa=0;aa<20;aa++) {
			ResBBRotamer* rot = new ResBBRotamer(s, atLib);
			rot->setAAType(aa, atLib);
			this->allRotLib1w[aa][id] = rot;
		}
		id++;
	}
	file.close();


	libFile = path + "/bbRotamer/rot1w/all.index";
	file.open(libFile.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open rotamer index file " << libFile << endl;
		exit(1);
	}
	id = 0;
	while(getline(file, s)){
		this->index1wTo1k[id] = atoi(s.c_str());
		id++;
	}
	file.close();

}

ResBBRotamerLib::~ResBBRotamerLib() {
	// TODO Auto-generated destructor stub


	delete this->atLib;

	for(int aa=0;aa<20;aa++) {
		for(int i=0;i<1000;i++){
			delete allRotLib1k[aa][i];
		}
	}

	for(int aa=0;aa<20;aa++) {
		for(int i=0;i<10000;i++){
			delete allRotLib1w[aa][i];
		}
	}


}

}
