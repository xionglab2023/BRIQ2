/*
 * RNARMSD.cpp
 *
 *  Created on: 2022��5��10��
 *      Author: pengx
 */

#include "geometry/RMSD.h"
#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>

#include "model/StructureModel.h"

using namespace NSPmodel;
using namespace std;


void printHelp(){
	cout << "Usage: rna_rms $RNAPDB1 $RNAPDB2" << endl;
}

int main(int argc, char** argv){

	if(argc != 3 || argv[1] == "-h"){
		printHelp();
		exit(0);
	}

	string pdbfile1 = string(argv[1]);
	string pdbfile2 = string(argv[2]);

	RNAPDB *pdb1 = new RNAPDB(pdbfile1, "");
	RNAPDB *pdb2 = new RNAPDB(pdbfile2, "");

	vector<RNABase*> baseList1 = pdb1->getBaseList();
	vector<RNABase*> baseList2 = pdb2->getBaseList();

	int n1 = baseList1.size();
	int n2 = baseList2.size();
	if(n1 != n2){
		cout << "sequence length not equal" << endl;
		exit(1);
	}

	vector<XYZ> tList1;
	vector<XYZ> tList2;
	for(int i=0;i<n1;i++){
		Atom* a1 = baseList1[i]->getAtom("C1'");
		Atom* a2 = baseList2[i]->getAtom("C1'");
		if(a1 == NULL || a2 == NULL) {
			cout << "base " << i <<  " miss C1' atom: " << endl;
			exit(0);
		}
		tList1.push_back(baseList1[i]->getAtom("C1'")->coord);
		tList2.push_back(baseList2[i]->getAtom("C1'")->coord);
	}

	double rms = rmsd(tList1, tList2);
	printf("%6.4f\n", rms);

}
