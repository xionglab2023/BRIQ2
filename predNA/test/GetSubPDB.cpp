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


void getSubPDB(string& pdbFile, string& subSeq, string outFile){
    RNAPDB* pdb = new RNAPDB(pdbFile, "xxxx");
    AtomLib* atLib = new AtomLib();
    vector<RNABase*> baseList = pdb->getValidBaseList(atLib);
    if(baseList.size() != subSeq.length()) {
        cout << "base number is not equal to subSeq length: " << baseList.size() << " " << subSeq.length() << endl;
        exit(0);
    }
    RNAChain* rc = new RNAChain("A");
    for(int i=0;i<baseList.size();i++){
        if(subSeq[i] == '-') continue;
        rc->addBase(baseList[i]);
    }
    ofstream out;
    out.open(outFile.c_str(), ios::out);
    if(!out.is_open()) 
    {
        cout << "fail to open file: " << outFile << endl;
    }
    rc->printPDBFormat(out, 1);
    out.close();

    delete rc;
    delete pdb;
    delete atLib;
}

int main(int argc, char** argv){

	if(argc !=  4 || argv[1] == "-h"){
		cout << "Usage: subPDB $PDBFILE $subSeq $OUTPDBFILE" << endl;
		exit(0);
	}
    string pdbFile = string(argv[1]);
    string subSeq = string(argv[2]);
    string outFile = string(argv[3]);
    getSubPDB(pdbFile, subSeq, outFile);

}