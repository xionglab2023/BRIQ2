/*
 * CalculateDetailEnergy.cpp
 *
 *  Created on: 2022Äê1ÔÂ29ÈÕ
 *      Author: pengx
 */

#include "protein/SeqDesignTemplateFast.h"
#include "para/ProParameter.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>

using namespace std;
using namespace NSPprotein;
using namespace NSPpara;
using namespace NSPtools;

int main(int argc, char** argv){

	//string filePolar = "/user/xiongpeng/cpp/ProteinModeling/data/para/pDes";
	//string designPara = "/user/xiongpeng/cpp/ProteinModeling/data/para/design/pz";

	string pdbID = string(argv[1]);
	string rankFile = string(argv[2]);

	clock_t start = clock();
	srand(time(0));
	cout << "init pa" << endl;
	ProParameter* pa = new ProParameter();
	DesignPara* dp = new DesignPara("/user/xiongpeng/sword/designPara/bp4/p10");
	cout << "init et" << endl;
	ProEnergyTable* pe = new ProEnergyTable();
	cout << "init s1s2" << endl;
	ProS1S2Energy* eS1S2= new ProS1S2Energy(dp, 0);
	cout << "init ec" << endl;
	EnergyCalculator* ec = new EnergyCalculator(pe, dp, pa);

	string inputFile = "/user/xiongpeng/sword/train/"+pdbID+"/input";
	string s2File = "/user/xiongpeng/sword/train/S2/"+pdbID+".s2";
	ResBBRotamerLib* bbLib = new ResBBRotamerLib();
	ResScRotamerLibMini* scLib = new ResScRotamerLibMini();
	AtomLib* atLib = new AtomLib();

	cout << "init" << endl;
	SeqDesignTemplateFast* dt = new SeqDesignTemplateFast(inputFile, bbLib, scLib, atLib, ec, eS1S2, s2File);
	dt->printDetailEnergy();

	ofstream out;
	out.open(rankFile.c_str(), ios::out);
	dt->getPositionRanks(out);
	out.close();


}


