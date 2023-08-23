/*
 * CalRNAEnergy.cpp
 *
 */




#include "geometry/localframe.h"
#include "geometry/TransMatrix.h"
#include "model/BaseDistanceMatrix.h"
#include "model/BasePairComposition.h"
#include <time.h>
#include <stdio.h>

#include "model/StructureModel.h"
#include "forcefield/RnaEnergyTable.h"

using namespace NSPmodel;
using namespace NSPforcefield;


int main(int argc, char** argv){

	clock_t start = clock();
	string pdbFile = string(argv[1]);
	RNAPDB pdb = RNAPDB(pdbFile, "xxxx");
	vector<RNABase*> baseList = pdb.getBaseList();
	RnaEnergyTable et;
	double tot = 0;
	printf("total energy: %7.3f\n",tot);

}
