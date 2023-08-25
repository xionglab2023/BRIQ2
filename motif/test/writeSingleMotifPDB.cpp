/*
 * writeSingleMotifPDB.cpp
 *
 *  Created on: 2023Äê8ÔÂ25ÈÕ
 *      Author: nuc
 */


#include "motif/RNAGraph.h"
#include "model/StructureModel.h"
#include <stdlib.h>
#include <iostream>
#include <time.h>

using namespace NSPmotif;
using namespace std;
using namespace NSPmodel;

int main(int argc, char** argv){

	string pdbFile = string(argv[1]);
	string mtfSeq = string(argv[2]);
	string outFile = string(argv[3]);
	RNAPDB* pdb = new RNAPDB(pdbFile);
	AtomLib* atLib = new AtomLib();
	vector<RNABase*> baseList = pdb->getValidBaseList(atLib);
	RNAPDB* rp = new RNAPDB();
	RNAChain* rc = NULL;
	string lastChainID = "unk";
	for(int k=0;k<mtfSeq.length();k++){
		string curChainID = mtfSeq.substr(k,1);
		if(mtfSeq[k] == '.') continue;
		if(curChainID == lastChainID){
			baseList[k]->chainID = curChainID;
			rc->addBase(baseList[k]);
		}
		else {
			baseList[k]->chainID = curChainID;
			lastChainID = curChainID;
			rc = new RNAChain();
			rc->setChainID(lastChainID);

			rp->addChain(rc);
			rc->addBase(baseList[k]);
		}
	}
	int baseNum = rp->getBaseList().size();
	ofstream out;
	out.open(outFile.c_str(), ios::out);
	rp->printPDBFormat(out);
	out.close();
}


