/*
 * TestBaseAtomNameAndSep.cpp
 *
 *  Created on: 2023Äê6ÔÂ8ÈÕ
 *      Author: nuc
 */


#include "geometry/localframe.h"
#include "geometry/TransMatrix.h"
#include "model/BaseDistanceMatrix.h"
#include "model/BasePairComposition.h"
#include "forcefield/RiboseOxygenEnergyTable.h"
#include "forcefield/XPara.h"
#include "model/BaseRotamer.h"
#include "model/RiboseRotamer.h"
#include <time.h>
#include <stdio.h>
#include <vector>
#include "forcefield/ForceFieldPara.h"
#include "forcefield/AtomicClashEnergy.h"
#include "model/StructureModel.h"
#include "forcefield/RnaAtomicEnergyTable.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace std;

int main(int argc, char** argv){

	ForceFieldPara* ffp = new ForceFieldPara();
	AtomicClashEnergy* ace = new AtomicClashEnergy(ffp);

	AtomLib* atLib = new AtomLib();

	vector<string> types;
	types.push_back("A");
	types.push_back("U");
	types.push_back("G");
	types.push_back("C");
	types.push_back("DA");
	types.push_back("DT");
	types.push_back("DG");
	types.push_back("DC");

	for(int k=0;k<8;k++){
		vector<string>* names = atLib->getRnaAtomNames(k);
		vector<string>* baseNames = atLib->getRnaSidechainAtoms(k);
		cout << endl;
		for(int i=0;i<baseNames->size();i++){
			string atomName = baseNames->at(i);
			string uniqueName = types[k]+"-"+atomName;
			int uniqueID = atLib->uniqueNameToUniqueID[uniqueName];
			printf("%-7s %2d %3d\n", uniqueName.c_str(), i, uniqueID);
			//cout << baseNames->at(i) << endl;
		}
	}




}

