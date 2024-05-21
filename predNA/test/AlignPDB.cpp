#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>

#include "model/StructureModel.h"
#include "predNA/BRNode.h"
#include "predNA/BRFoldingTree.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredna;
using namespace std;


void align(string alignedTarget,  string inputFile, string outFile){

    vector<XYZ> tList1;
	vector<XYZ> tList2;



    RNAPDB* pdbA = new RNAPDB(alignedTarget);
    RNAPDB* pdbB = new RNAPDB(inputFile);

    vector<RNABase*> baseListA = pdbA->getBaseList();
    vector<RNABase*> baseListB = pdbB->getBaseList();

    int seqLen = baseListA.size();

    if(baseListA.size() != baseListB.size()){
        cout << "pdb length not equal: " << baseListA.size() << " " << baseListB.size() << endl;
    }

    for(int i=0;i<baseListA.size();i++){
        
        tList1.push_back(baseListA[i]->getAtom("C1'")->getCoord());
        tList2.push_back(baseListB[i]->getAtom("C1'")->getCoord());

    }

    double rms = rmsd(tList1, tList2);
    cout << "rms: " << rms << endl;

	XYZ Acog = getCOG(tList1);
	XYZ Bcog = getCOG(tList2);

    cout << Acog.toString() << endl;
    cout << Bcog.toString() << endl;


	vector<XYZ> listA;
	vector<XYZ> listB;
	for(int i=0;i<tList1.size();i++){
		XYZ a = tList1[i] - Acog;
		XYZ b = tList2[i] - Bcog;
		listA.push_back(a);
		listB.push_back(b);
	}
	TransForm tf = buildRotation(listA, listB);

    tf.print();

	string s = "AUGC";
	char ss[20];

	for(int i=0;i<seqLen;i++) {

		RNABase* base = baseListB[i];
		vector<Atom*>* aList = base->getAtomList();
		for(int j=0;j<aList->size();j++){
			Atom* a = aList->at(j);
            XYZ t = a->coord - Bcog;
			XYZ newCoord = tf.transform(t);
			newCoord = newCoord + Acog;
        
			a->setCoord(newCoord);
		}
	}



	ofstream of;
	of.open(outFile.c_str(), ios::out);
    pdbB->printPDBFormat(of);
	of.close();

	
}


int main(int argc, char** argv){

	if(argc !=  4 || argv[1] == "-h"){
		cout << "Usage: alignPDB $PDBFILE $targetFile $OUTPDBFILE" << endl;
		exit(0);
	}

    align(string(argv[1]), string(argv[2]), string(argv[3]));

}