/*
 * TestProPolarEnergy.cpp
 *
 */

#include <time.h>
#include <iostream>
#include "forcefield/ProAtomicEnergyTable.h"

using namespace std;
using namespace NSPforcefield;

int main(int argc, char** argv){

	cout << "testPE: " << endl;
	clock_t start = clock();
	cout << "init polar energy table: " << endl;

	//ProPolarEnergy* pe = new ProPolarEnergy();
	clock_t end = clock();
	cout << "time: " << (float)(end-start)/CLOCKS_PER_SEC << "s" << endl;

	string filePolar = "/user/xiongpeng/cpp/ProteinModeling/data/para/pDes";
	string designPara = "/user/xiongpeng/cpp/ProteinModeling/data/para/design/px";

	ProParameter* pa = new ProParameter(1, designPara, filePolar);

	ProAtomicEnergyTable* at = new ProAtomicEnergyTable(pa);
	for(double d = 0.1;d<8.0;d+=0.02)
	{
		double d0 = 3.9;
		double shift = 0.0;
		double wd = -1.0;
		double lamda = 2.4;
		double e = at->vdwEnergyDesign(d, d0, shift, wd, lamda, 6.5);
		printf("%4.2f %8.3f\n", d, e);
	}

}



