/*
 * RiboseRotamerLib.cpp
 *
 *  Created on: 2022Äê9ÔÂ6ÈÕ
 *      Author: pengx
 */

#include "model/RiboseRotamerLib.h"

namespace NSPmodel {


RiboseRotamerLib::RiboseRotamerLib() {
	// TODO Auto-generated constructor stub
	string augc = "AUGCatgc";
	string path = NSPdataio::datapath();

	ifstream file;
	string s;
	XPara para;
	for(int i=0;i<8;i++){
		string baseType = augc.substr(i,1);
		string fileName = path+"riboseRotamer/"+baseType+".rot";
		cout << "init " << fileName << endl;
		file.open(fileName.c_str(), ios::in);
		if(! file.is_open()) {
			cout << "fail to open rotamer library file " << fileName << endl;
			exit(1);
		}
		int index = 0;
		while(getline(file,s)){
			RiboseRotamer* rot = new RiboseRotamer(s);

			rot->resType = i;
			rot->rotTypeLv1 = 0;
			rot->rotType = index;
			rotLib[i][index] = rot;
			index++;
		}
		file.close();
	}
}

RiboseRotamerLib::RiboseRotamerLib(XPara* para){
	string augc = "AUGCatgc";
	string path = NSPdataio::datapath();
	ifstream file;
	string s;
	for(int i=0;i<8;i++){
		string baseType = augc.substr(i,1);
		string fileName = path+"riboseRotamer/"+baseType+".rot";
		file.open(fileName.c_str(), ios::in);
		if(! file.is_open()) {
			cout << "fail to open rotamer library file " << fileName << endl;
			exit(1);
		}

		int index = 0;
		while(getline(file,s)){
			RiboseRotamer* rot = new RiboseRotamer(s);
			rot->energy = rot->energy * para->wtRotamer;
			rot->resType = i;
			rot->rotTypeLv1 = 0;
			rot->rotType = index;
			rotLib[i][index] = rot;
			index++;
		}
		file.close();

	}
}


RiboseRotamerLib::~RiboseRotamerLib() {
	// TODO Auto-generated destructor stub
	for(int i=0;i<8;i++){
		for(int j=0;j<1500;j++){
			delete rotLib[i][j];
		}
	}
}

} /* namespace NSPmodel */
