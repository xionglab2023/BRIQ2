/*
 * TestRotLib.cpp
 *
 *  Created on: 2023Äê12ÔÂ1ÈÕ
 *      Author: nuc
 */



#include "model/RNABaseLib.h"
#include "geometry/RMSD.h"
#include "stdlib.h"
#include <iostream>
#include "model/StructureModel.h"
#include "model/RiboseRotamerLib.h"

using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;

int main(int argc, char** argv){

	RNAPDB* pdb = new RNAPDB(string(argv[1]), "xxxx");


	vector<RNABase*> baseList = pdb->getBaseList();
	cout << baseList.size()	<< endl;

	RNABase* base = baseList[3];

	RiboseRotamerLib* rotLib = new RiboseRotamerLib();
	RiboseRotamer* rot = rotLib->getNearestRotamer(base);
	cout << "chi: " << rot->chi << endl;
	cout << "imp: " << rot->improper << endl;
	cout << "ene: " << rot->energy << endl;
	delete pdb;
	delete rotLib;


}



