/*
 * BuildMutation.cpp
 *
 *  Created on: 2022Äê5ÔÂ11ÈÕ
 *      Author: pengx
 */

#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>

#include "model/StructureModel.h"
#include "predNA/BRNode.h"
#include "predNA/BRFoldingTree.h"
#include "predNA/MCRun.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredna;
using namespace std;


int main(int argc, char** argv){

	if(argc !=  4 || argv[1] == "-h"){
		cout << "Usage: briq_buildMutation $PDBFILE $SEQUENCE $OUTPDBFILE" << endl;
		exit(0);
	}

	RotamerLib* rotLib = new RotamerLib();

	string pdbfile = string(argv[1]);
	string seq = string(argv[2]);
	string outpdb = string(argv[3]);
	int len = seq.length();

	RNAPDB pdb(pdbfile, "xxx");
	vector<RNABase*> baseList= pdb.getBaseList();
	if(baseList.size() != seq.length()) {
		cout << "sequence length not equal to pdb file" << endl;
		exit(1);
	}

	map<char, int> c2i;
	c2i['A'] = 0;
	c2i['U'] = 1;
	c2i['G'] = 2;
	c2i['C'] = 3;
	int typeList[len];
	vector<LocalFrame> csList;
	bool connectToNeighbor[len];

	for(int i=0;i<seq.length();i++){
		char c = seq[i];
		if(c != 'A' && c!= 'U' && c!= 'G' && c!= 'C') {
			cout << "invalid sequence "<< seq << endl;
			exit(0);
		}

		baseList[i]->updateCoordSystem();
		if(!baseList[i]->hasLocalFrame) {
			cout << "base " << i << " atom not complete" << endl;
			exit(0);
		}
		LocalFrame cs = baseList[i]->getCoordSystem();
		typeList[i] = c2i[c];
		csList.push_back(cs);
	}
	for(int i=0;i<len-1;i++){
		if(baseList[i]->connectToNeighbor(baseList[i+1]))
			connectToNeighbor[i] = true;
		else
			connectToNeighbor[i] = false;
	}
	connectToNeighbor[len-1] = false;

	BRNode** nodeList = new BRNode*[len];

	//
	cout << "generate nodes" << endl;
	for(int i=0;i<len;i++){
		cout << "node: " << i << endl;
		BRNode* node = new BRNode(typeList[i], i, rotLib);
		node->baseConf->updateCoords(csList[i]);

		RiboseRotamer* rot;

		if(!baseList[i]->backboneComplete())
			rot = rotLib->riboseRotLib->getLowestEnergyRotamer(baseList[i]->baseTypeInt);
		else
			rot = new RiboseRotamer(baseList[i]);

		node->riboseConf->updateLocalFrame(csList[i]);
		node->riboseConf->updateRotamer(rot);

		nodeList[i] = node;
	}

	cout << "build po3" << endl;

	ForceFieldPara* para = new ForceFieldPara();
	PO3Builder* pb = new PO3Builder(para);
	for(int i=0;i<len-1;i++){
		BRNode* nodeA = nodeList[i];
		BRNode* nodeB = nodeList[i+1];
		nodeA->connectToNeighbor = true;
		pb->buildPhosphate(nodeA->riboseConf, nodeB->riboseConf, nodeA->phoConf);
	}

	BRTreeInfo* info = new BRTreeInfo(len, typeList, connectToNeighbor, nodeList, 0.0, rotLib);
	info->printPDB(outpdb);

}

