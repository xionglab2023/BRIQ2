/*
 * GenerateInitModel.cpp
 *
 */

#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>

#include "model/StructureModel.h"
#include "predNA/BRFoldingTree.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredNA;
using namespace std;


int main(int argc, char** argv){

	if(argc !=  3 || argv[1] == "-h"){
		cout << "Usage: briq_initPDB $SEQUENCE $OUTPDBFILE" << endl;
		exit(0);
	}
	RotamerLib* rotLib = new RotamerLib();

	string seq = string(argv[1]);
	string outpdb = string(argv[2]);

	CsMove cm("1.179875   3.813016   -3.749098  0.8472664  0.5174403  -0.1199800 -0.5178250 0.8549461  0.0304035  0.1183084  0.0363688  0.9923107");

	vector<int> typeList;
	vector<LocalFrame> csList;

	map<char, int> c2i;
	c2i['A'] = 0;
	c2i['U'] = 1;
	c2i['G'] = 2;
	c2i['C'] = 3;

	LocalFrame cs;

	csList.push_back(cs);
	typeList.push_back(c2i[seq[0]]);

	for(int i=1;i<seq.length();i++){
		cs = cs + cm;

		char c = seq[i];
		if(c2i.find(c) == c2i.end())
			continue;
		else {
			cout << c << endl;
			typeList.push_back(c2i[c]);
			csList.push_back(cs);
		}
	}

	cout << "a" << endl;

	int len = typeList.size();
	bool connectToNeighbor[len];
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

	cout << "b" << endl;

	id = 0;
	for(int i=1;i<seq.length();i++){
		if(c2i.find(seq[i]) == c2i.end()){
			connectToNeighbor[id] = false;
			id++;
			i++;
		}
		else {
			connectToNeighbor[id] = true;
			id++;
		}
	}
	connectToNeighbor[len-1] = false;

	cout << "c" << endl;

	BRNode** nodeList = new BRNode*[len];
	for(int i=0;i<len;i++){
		BRNode* node = new BRNode(typeList[i], i, rotLib);
		node->baseConf->updateCoords(csList[i]);
		node->baseConfTmp->updateCoords(csList[i]);
		node->riboseConf->updateLocalFrameAndRotamer(csList[i], rotLib->riboseRotLib->getLowestEnergyRotamer(0));
		node->riboseConf->updateLocalFrame(csList[i]);
		node->riboseConfTmp->copyValueFrom(node->riboseConf);
		nodeList[i] = node;
	}

	ForceFieldPara* para = new ForceFieldPara();
	PO3Builder* pb = new PO3Builder(para);
	for(int i=0;i<len-1;i++){

		BRNode* nodeA = nodeList[i];
		BRNode* nodeB = nodeList[i+1];
		nodeA->connectToNeighbor = connectToNeighbor[i];
		pb->buildPhosphate(nodeA->riboseConf, nodeB->riboseConf, nodeA->phoConf);
	}

	cout << "e " << endl;

	BRTreeInfo* info = new BRTreeInfo(len, baseTypeSeq, connectToNeighbor, nodeList, 0.0, rotLib);
	info->printPDB(outpdb);

	cout << "f " << endl;
}


