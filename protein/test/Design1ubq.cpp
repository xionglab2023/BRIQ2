/*
 * Design1ubq.cpp
 *
 *  Created on: 2022Äê7ÔÂ26ÈÕ
 *      Author: pengx
 */
#include "protein/SeqDesignTemplate.h"
#include "protein/SeqDesignTemplateFast.h"
#include "para/ProParameter.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>

using namespace std;
using namespace NSPprotein;
using namespace NSPpara;
using namespace NSPtools;

void design1ubq(){
	clock_t start = clock();
	srand(time(0));
	ProParameter* pa = new ProParameter();
	string paraFile = "/user/xiongpeng/sword/designPara/bp/p20";


	cout << paraFile << endl;
	DesignPara* dp = new DesignPara(paraFile);
	ProEnergyTable* pe = new ProEnergyTable();
	ProS1S2Energy* eS1S2= new ProS1S2Energy(dp, 0);



	EnergyCalculator* ec = new EnergyCalculator(pe, dp, pa);

	string inputFile = "/user/xiongpeng/input";
	string output = "/user/xiongpeng/1ubq-desSeq";
	string outpdb = "/user/xiongpeng/1ubq.des.pdb";
	ofstream of;
	of.open(output.c_str(), ios::out);

	ResBBRotamerLib* bbLib = new ResBBRotamerLib();
	ResScRotamerLibMini* scLib = new ResScRotamerLibMini();
	AtomLib* atLib = new AtomLib();
	SeqDesignTemplateFast* dt = new SeqDesignTemplateFast(inputFile, bbLib, scLib, atLib, ec, eS1S2);
	dt->loadEAEM();

	clock_t t1 = clock();
	cout << "time: " << (float)(t1-start)/CLOCKS_PER_SEC << "s" << endl;

	dt->printNatEnergy();


	dt->designMC();
	clock_t t2 = clock();
	cout << "time: "<< "des " << (float)(t2-start)/CLOCKS_PER_SEC << "s" << endl;
	of << dt->getDesignSequence() << endl;
	dt->printDesignPDB(outpdb);
	of.close();

	delete bbLib;
	delete scLib;
	delete atLib;
	delete dt;
}

int main(int argc, char** argv){
	design1ubq();
}
