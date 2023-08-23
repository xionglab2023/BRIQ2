/*
 * TestAtomLib.cpp
 *
 *  Created on: 2022Äê1ÔÂ28ÈÕ
 *      Author: pengx
 */

#include "model/AtomLib.h"
#include "model/BaseRotamerLib.h"

using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;

int main(){
	cout << "test atLib" << endl;
	cout << "init" << endl;
	AtomLib* atLib = new AtomLib();
	for(int i=0;i<334;i++){
		AtomProperty* ap = atLib->apList[i];
		double meanRadii = ap->vdwRadius;
		string name = ap->atomUniqueName;
		//printf("%-7s %5.3f\n", name.c_str(), meanRadii);
	}

	BaseRotamerLib* baseLib = new BaseRotamerLib(atLib);
	baseLib->testRotamerLib(atLib);

}




