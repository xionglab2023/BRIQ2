/*
 * BuildMutation.cpp
 *
 *  Created on: 2022��5��11��
 *      Author: pengx
 */

#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>

#include "model/StructureModel.h"
#include "predNA/NuGraph.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredNA;
using namespace std;


int main(int argc, char** argv){

	if(argc !=  4 || argv[1] == "-h"){
		cout << "Usage: rna_buildMutation $PDBFILE $SEQUENCE $OUTPDBFILE" << endl;
		exit(0);
	}

	RotamerLib* rotLib = new RotamerLib();
	AtomLib* atLib = new AtomLib();
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
	bool fixed[len];
	

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
		fixed[i] = false;
		if(baseList[i]->connectToNeighbor(baseList[i+1]))
			connectToNeighbor[i] = true;
		else
			connectToNeighbor[i] = false;
	}
	connectToNeighbor[len-1] = false;

	NuNode** nodeList = new NuNode*[len];

	//
	cout << "generate nodes" << endl;
	for(int i=0;i<len;i++){
		cout << "node: " << i << endl;
		BaseRotamer* baseRot = rotLib->baseRotLib->baseLib[typeList[i]];
		

		RiboseRotamer* rot;

		if(!baseList[i]->backboneComplete())
			rot = rotLib->riboseRotLib->getLowestEnergyRotamer(baseList[i]->baseTypeInt);
		else
			rot = new RiboseRotamer(baseList[i]);

		NuNode* node = new NuNode(i, typeList[i], csList[i], baseRot, rot, atLib);
		node->riboseConf->updateLocalFrameAndRotamer(csList[i], rot);
		nodeList[i] = node;
	}

	cout << "build po3" << endl;

	ForceFieldPara* para = new ForceFieldPara();
	PO3Builder* pb = new PO3Builder(para);
	for(int i=0;i<len-1;i++){
		NuNode* nodeA = nodeList[i];
		NuNode* nodeB = nodeList[i+1];
		nodeA->connectToNeighbor = true;
		pb->buildPhosphate(nodeA->riboseConf, nodeB->riboseConf, nodeA->phoConf);
	}

	//graphInfo* info = new graphInfo(len, )
	graphInfo* info = new graphInfo(len, typeList, connectToNeighbor, fixed, nodeList, 0.0, atLib, 0);
	info->printPDB(outpdb);

}

