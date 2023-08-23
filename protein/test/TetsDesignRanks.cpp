/*
 * TestDesign.cpp
 *
 *  Created on: 2022Äê1ÔÂ9ÈÕ
 *      Author: pengx
 */

#include "protein/SeqDesignTemplate.h"
#include "para/ProParameter.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>

using namespace std;
using namespace NSPprotein;
using namespace NSPpara;
using namespace NSPtools;

void test1rr2(){
	clock_t start = clock();
	srand(time(0));
	ResName rn;

	ProParameter* pa = new ProParameter();
	DesignPara* dp = new DesignPara("/user/xiongpeng/sword/designPara/py");
	ProEnergyTable* pe = new ProEnergyTable();
	ProS1S2Energy* eS1S2= new ProS1S2Energy(dp, 0);

	EnergyCalculator* ec = new EnergyCalculator(pe, dp, pa);
	string inputFile = "/user/xiongpeng/input.1rr2";
	ResBBRotamerLib* bbLib = new ResBBRotamerLib();
	ResScRotamerLib* scLib = new ResScRotamerLib();
	AtomLib* atLib = new AtomLib();
	SeqDesignTemplate* dt = new SeqDesignTemplate(inputFile, bbLib, scLib, atLib, ec, eS1S2);
	cout << "get ranks" << endl;
	//dt->printS1S2();
	//dt->printPairInfo();
	//dt->getPositionRanks();

}

int main(int argc, char** argv){

	test1rr2();

}
