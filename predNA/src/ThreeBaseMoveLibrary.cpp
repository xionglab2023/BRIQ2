/*
 * ThreeBaseMoveLibrary.cpp
 *
 */

#include "predNA/ThreeBaseMoveLibrary.h"

namespace NSPpredna {

ThreeBaseMoveLibrary::ThreeBaseMoveLibrary(int typeA, int typeB, int typeC) {
	// TODO Auto-generated constructor stub
	string path = NSPdataio::datapath()+"rna/";
	string libPath = path + "moveLib/threeBase";
	ifstream file;
	this->rotNum1 = 50;
	this->rotNum2 = 100;
	string augc = "AUGC";

	string tag = ""+augc.substr(typeA,1)+augc.substr(typeB,1)+augc.substr(typeC,1);
	string fileA = libPath+"/acConnect/"+tag+".rot";
	string fileB = libPath+"/all/"+tag+".rot";
	file.open(fileA, ios::in);
	string s;
	getline(file, s);
	while(getline(file,s)){
		string s1 = s.substr(0,119);
		string s2 = s.substr(120, 119);
		string s3 = s.substr(240, 119);
		CsMove mv1(s1);
		CsMove mv2(s2);
		CsMove mv3(s3);
		mvLib1AB.push_back(mv1);
		mvLib1BC.push_back(mv2);
	}
	file.close();

	file.open(fileB, ios::in);

	getline(file, s);
	while(getline(file,s)){
		string s1 = s.substr(0,119);
		string s2 = s.substr(120, 119);
		string s3 = s.substr(240, 119);
		CsMove mv1(s1);
		CsMove mv2(s2);
		CsMove mv3(s3);
		mvLib2AB.push_back(mv1);
		mvLib2BC.push_back(mv2);
	}
	file.close();
}

ThreeBaseMoveLibrary::~ThreeBaseMoveLibrary() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPforcefield */
