/*
 * SinglePackToPDB.cpp
 *
 */


#include "protein/SingleSidechainModelingTemplate.h"
#include "para/ProParameter.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>

using namespace std;
using namespace NSPprotein;
using namespace NSPpara;
using namespace NSPtools;

void writePDB(){


	ResName* rn = new ResName();
	ResScRotamerLib* scLib = new ResScRotamerLib();
	ResBBRotamerLib* bbLib = new ResBBRotamerLib();
	AtomLib* atLib = new AtomLib();
	EnergyCalculator* ec = NULL;
	string s;
	vector<string> spt;
	vector<string> lines;

	for(int aa=0;aa<20;aa++){
		string tri = rn->intToTri(aa);
		if(tri == "GLY") continue;
		string datFile = "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/split/" + tri + ".dat";
		ifstream file;
		file.open(datFile, ios::in);
		if(!file.is_open()){
			cout << "fail to open file: " << datFile << endl;
			exit(1);
		}
		int index = 1;
		char xx[200];
		while(getline(file, s)){
			if(s[0] == '#'){
				splitString(s, " ", &spt);
				int neighborNum = atoi(spt[2].c_str());
				lines.clear();
				for(int i=0;i<neighborNum+1;i++){
					getline(file, s);
					lines.push_back(s);
				}
				SingleSidechainModelingTemplate* st = new SingleSidechainModelingTemplate(lines, bbLib, scLib, rn);


				sprintf(xx, "%d", index);
				string sid = string(xx);
				string out1 = "/user/xiongpeng/proteinModeling/singleSidechainModeling/pdb/" + tri + "-" + sid + ".pdb";
				st->printInit(out1, atLib);
				index++;

				delete st;
				if(index == 20) break;
			}
		}
		file.close();
	}

	delete rn;
	delete scLib;
	delete bbLib;
}


int main(int argc, char** argv){

	writePDB();
}

