/*
 * TestMerge.cpp
 *
 */

#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include "forcefield/XPara.h"
#include "model/StructureModel.h"
#include "predNA/BRFoldingTree.h"
#include "predNA/MCRun.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredna;
using namespace std;

int main(int argc, char** argv){

	RiboseRotamerLib rotLib;

	string file1 = string(argv[1]);
	string file2 = string(argv[2]);
	string outpdb = string(argv[3]);

	RNAPDB* frag1 = new RNAPDB(file1, "xxxx");
	RNAPDB* frag2 = new RNAPDB(file2, "xxxx");

	vector<RNABase*> baseListA = frag1->getBaseList();
	vector<RNABase*> baseListB = frag2->getBaseList();

	CsMove cm1("1.179875   3.813016   -3.749098  0.8472664  0.5174403  -0.1199800 -0.5178250 0.8549461  0.0304035  0.1183084  0.0363688  0.9923107");
	LocalFrame csA;
	LocalFrame csB;
	csB = csA + cm1;
	CsMove cm2 = csA - csB;

	LocalFrame csList[39];

	csList[2] = baseListA[0]->getCoordSystem() + cm2;
	csList[1] = csList[2]+ cm2;
	csList[0] = csList[1] + cm2;

	csList[36] = baseListA[20]->getCoordSystem() + cm1;
	csList[37] = csList[36] + cm1;
	csList[38] = csList[37] + cm1;

	for(int i=0;i<11;i++){
		csList[i+3] = baseListA[i]->getCoordSystem();
	}

	for(int i=0;i<10;i++){
		csList[i+26] = baseListA[i+11]->getCoordSystem();
	}
	csList[14] = csList[13] + cm1;
	csList[15] = csList[14] + cm1;
	csList[16] = csList[15] + cm1;



	vector<LocalFrame> frag2CsList;
	for(int i=0;i<baseListB.size();i++){
		frag2CsList.push_back(baseListB[i]->getCoordSystem());
	}

	for(int i=0;i<7;i++){
		CsMove mv = frag2CsList[i+1] - frag2CsList[0];
		csList[17+i] = csList[16] + mv;
	}

	csList[23].print();

	csList[24] = csList[23] + cm1;

	csList[24].print();

	csList[25] = csList[24] + cm1;
	csList[25].print();

	string seq = "GGCAGAGCUCAACACAGCGAAAGCUGUGGCUAGACUGUC";
	map<char, int> c2i;
	c2i['A'] = 0;
	c2i['U'] = 1;
	c2i['G'] = 2;
	c2i['C'] = 3;

	int len = seq.length();
	int baseTypeSeq[len];
	int id = 0;
	for(int i=0;i<seq.length();i++){
		if(c2i.find(seq[i]) == c2i.end()){
			continue;
		}
		else{
			baseTypeSeq[id] = c2i[seq[i]];
			id++;
		}
	}
	bool connectToNeighbor[len];
	for(int i=0;i<len;i++){
		connectToNeighbor[i] = true;
	}

	/*
	BRNode** nodeList = new BRNode*[len];
	for(int i=0;i<len;i++){
		BRNode* node = new BRNode(baseTypeSeq[i], i);
		node->cs1 = csList[i];
		node->rot = rotLib.rotLib[2][163];
		node->cs2 = node->cs1 + node->rot->mv12;
		for(int j=0;j<node->baseAtomNum;j++){
			node->baseAtomCoords[j] = local2global(node->cs1, node->atomCoordLocal[j]);
			node->baseAtomCoordsTmp[j] = node->baseAtomCoords[j];
		}
		for(int j=0;j<11;j++){
			node->riboAtomCoords[j] = local2global(node->cs1, node->rot->tList1[j]);
			node->riboAtomCoordsTmp[j] = node->riboAtomCoords[j];
		}
		nodeList[i] = node;
	}

	XPara para;
	RnaAtomicEnergyTable* at = new RnaAtomicEnergyTable(&para);
	RiboseOxygenEnergyTable* rET = new RiboseOxygenEnergyTable();
	PO3Builder* pb = new PO3Builder(&para);
	for(int i=0;i<len-1;i++){

		BRNode* nodeA = nodeList[i];
		BRNode* nodeB = nodeList[i+1];
		nodeA->connectToNeighbor = true;
		nodeA->phoLocal = pb->getPhoLocal(nodeA->cs2, nodeB->riboAtomCoords, nodeA->rot->improper, nodeB->rot->improper);
		nodeA->pho = PhophateGroup(nodeA->phoLocal, nodeA->cs2);
	}

	vector<int> freeIDs;
	for(int i=0;i<len;i++){
		freeIDs.push_back(i);
	}
	*/


}


