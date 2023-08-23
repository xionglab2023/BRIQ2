/*
 * TestBaseLib.cpp
 *
 */


#include "model/RNABaseLib.h"
#include "geometry/RMSD.h"
#include "stdlib.h"
#include <iostream>
#include "model/StructureModel.h"


using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;

int main(int argc, char** argv){

	RNAPDB* pdb = new RNAPDB(string(argv[1]), "xxxx");
	vector<RNABase*> baseList = pdb->getBaseList();

	RNABase* base = baseList[3];
	LocalFrame cs = base->getCoordSystem();
	XYZ tn1 = base->getAtom("N9")->coord;
	XYZ local = global2local(cs, tn1);
	cout << local.toString() << endl;

}



