/*
 * TestAtomicEnergy.cpp
 *
 */

#include "protein/SingleSidechainModelingTemplate.h"
#include "para/ProParameter.h"
#include "model/ResScRotamerLib.h"
#include "model/ResBBRotamerLib.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include "model/StructureModel.h"

using namespace std;
using namespace NSPprotein;
using namespace NSPpara;
using namespace NSPtools;
using namespace NSPmodel;

int main(int argc, char** argv){

	ResName* rn = new ResName();
	ResScRotamerLib* scLib = new ResScRotamerLib();
	ResBBRotamerLib* bbLib = new ResBBRotamerLib();
	AtomLib* atLib = new AtomLib();

	ProParameter* pa1 = new ProParameter("/user/xiongpeng/proteinModeling/para/qSR");
	ProParameter* pa2 = new ProParameter("/user/xiongpeng/proteinModeling/para/qMR");
	ProParameter* pa3 = new ProParameter("/user/xiongpeng/proteinModeling/para/qLR");

	DesignPara* dp = new DesignPara("");
	ProEnergyTable* et = new ProEnergyTable();
	EnergyCalculator* ec = new EnergyCalculator(et, dp, pa1);

	string s;
	vector<string> spt;
	vector<string> lines;

	{
		string datFile = "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/arg.dat";
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
				cout << "init single sidechain modeling template" << endl;
				SingleSidechainModelingTemplate* st = new SingleSidechainModelingTemplate(lines, bbLib, scLib, rn);
				cout << "node1" << endl;
				ResNode* nd1 = st->targetNode;
				cout << "node2" << endl;
				ResNode* nd2 = st->neighborNodes[14];

				string type1 = rn->intToTri(nd1->aaType);
				string type2 = rn->intToTri(nd2->aaType);
				cout << rn->intToTri(nd1->aaType) << endl;
				cout << rn->intToTri(nd2->aaType) << endl;

				ResConformer* confA = nd1->conf;
				ResConformer* confB = nd2->conf;


				delete st;
			}
		}
		file.close();
	}

	delete rn;
	delete et;
	delete ec;
	delete scLib;
	delete bbLib;

}


