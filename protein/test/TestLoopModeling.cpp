/*
 * TestLoopModeling.cpp
 *
 */


#include "protein/LoopModelingTemplate.h"
#include "para/ProParameter.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>

using namespace std;
using namespace NSPprotein;
using namespace NSPpara;
using namespace NSPtools;

void printPDB(const string datFile){

	ifstream file;
	file.open(datFile, ios::in);
	if(!file.is_open()){
		cout << "fail to open file: " << datFile << endl;
		exit(1);
	}

	ResName* rn = new ResName();
	cout << "init scLib" << endl;
	ResScRotamerLib* scLib = new ResScRotamerLib();
	cout << "init bbLib" << endl;
	ResBBRotamerLib* bbLib = new ResBBRotamerLib();
	cout << "init atLib" << endl;
	AtomLib* atLib = new AtomLib();

	string s;
	vector<string> spt;
	vector<string> lines;
	vector<double> rmsList;
	int index = 1;
	char xx[200];
	while(getline(file, s)){
		if(s[0] == '#'){
			splitString(s, " ", &spt);
			int neighborNum = atoi(spt[1].c_str());
			lines.clear();
			for(int i=0;i<neighborNum+5;i++){
				getline(file, s);
				lines.push_back(s);
			}

			cout << "init LT " << index <<  endl;
			LoopModelingTemplate* lt = new LoopModelingTemplate(5, lines, bbLib, scLib, rn);

			sprintf(xx, "%d", index);

			string fileid = string(xx);
			index++;

			cout << "print pdb" << endl;
			string out = "/user/xiongpeng/proteinModeling/loopModeling/pdb/p" + fileid + "-init.pdb";
			lt->printInit(out, atLib);
			delete lt;
		}
	}

	delete rn;
	delete scLib;
	delete bbLib;
}

int main(int argc, char** argv){

	string dat = string(argv[1]);
	printPDB(dat);
}
