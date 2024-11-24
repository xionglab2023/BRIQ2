/*
 * ReplaceFragment.cpp
 *
 *  Created on: 2022��6��5��
 *      Author: pengx
 */

#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include "model/StructureModel.h"
#include "predNA/BRFoldingTree.h"
#include "forcefield/XPara.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredNA;
using namespace std;

int main(int argc, char** argv){

	if(argc !=  3 || argv[1] == "-h"){
		cout << "Usage: replaceFragment $INPUTFILE $OUTPUT" << endl;
		exit(0);
	}

	string inputFile = string(argv[1]);
	string output = string(argv[2]);

	RotamerLib* rotLib = new RotamerLib();

	NSPtools::InputParser input(inputFile);
	string pdbFile = input.getValue("pdb");
	string pdbFileM = input.getValue("pdbF");
	string baseSeq = input.getValue("seq");
	string fgA = input.getValue("fgA");

	int seqLen = baseSeq.length();
	BRNode** nodeList = new BRNode*[seqLen];

	RNAPDB pdb(pdbFile, "a");
	RNAPDB pdbF(pdbFileM, "b");
	vector<RNABase*> baseListP = pdb.getBaseList();
	vector<RNABase*> baseListF = pdbF.getBaseList();

	int fragLen = baseListF.size();

	cout << "parse fragment" << endl;

	if(fgA.length() != baseSeq.length()){
		cout << "fgA length error" << endl;
	}

	int nA = 0;

	int posX=-1;
	int posY=-1;
	map<int, int> seqPosToMtfPos;

	for(int i=0;i<baseSeq.length();i++){
		seqPosToMtfPos[i] = -1;
	}
	int m = 0;
	for(int i=0;i<baseSeq.length();i++){
		char c = fgA[i];
		if(c == '-')
			continue;

		seqPosToMtfPos[i] = m;
		m++;
	}
	if(m != fragLen) {
		cout << "fragment length not match: " << fragLen << " " << m << endl;
	}

	for(int i=0;i<baseSeq.length();i++){
		if(fgA[i] == 'A')
			nA++;
		else if(fgA[i] == 'X')
		{
			nA++;
			posX = i;
		}
		else if(fgA[i] != '-') {
			cout << "invalid fgA: " << fgA << endl;
			exit(0);
		}
	}

	if(posX < 0) {
		cout << "invalid fragment " << fgA << endl;
		exit(0);
	}


	map<char, int> c2i;
	c2i['A'] = 0;
	c2i['U'] = 1;
	c2i['G'] = 2;
	c2i['C'] = 3;
	int* typeList = new int[seqLen];
	LocalFrame csListP[seqLen];
	LocalFrame csListF[fragLen];

	LocalFrame csListA[seqLen];

	RNABase** baseList = new RNABase*[seqLen];

	cout << "read sequence" << endl;
	for(int i=0;i<seqLen;i++){
		char c = baseSeq[i];
		if(c != 'A' && c!= 'U' && c!= 'G' && c!= 'C') {
			cout << "invalid sequence "<< baseSeq << endl;
			exit(0);
		}
		typeList[i] = c2i[c];
	}

	int id = 0;
	for(int i=0;i<seqLen;i++){
		LocalFrame cs = baseListP[i]->getCoordSystem();
		baseList[i] = baseListP[id];
		csListP[i] = cs;
		csListA[i] = cs;
	}
	for(int i=0;i<fragLen;i++){
		LocalFrame cs = baseListF[i]->getCoordSystem();
		csListF[i] = cs;
	}

	for(int i=0;i<seqLen;i++){
		if(fgA[i] == '-' || i==posX){
			LocalFrame csOld = csListP[i];
			csListA[i] = csListP[i];
		}
		else {
			LocalFrame csOld = csListP[i];
			CsMove mvXF = csListF[seqPosToMtfPos[i]] - csListF[seqPosToMtfPos[posX]];
			LocalFrame csNew = csListP[posX] + mvXF;
			csListA[i] = csNew;
		}
	}

	cout << "generate nodes" << endl;

	for(int i=0;i<seqLen;i++){
		cout << "node: " << i << endl;
		BRNode* node = new BRNode(typeList[i], i, rotLib);
		node->baseConf->updateCoords(csListA[i]);

		RiboseRotamer* rot;
		if(!baseList[i]->backboneComplete())
			rot = rotLib->riboseRotLib->getLowestEnergyRotamer(baseList[i]->baseTypeInt);
		else
			rot = new RiboseRotamer(baseList[i]);

		node->riboseConf->updateLocalFrame(csListA[i]);
		node->riboseConf->updateRotamer(rot);
		nodeList[i] = node;
	}

	cout << "build po3" << endl;

	ForceFieldPara* para = new ForceFieldPara();
	PO3Builder* pb = new PO3Builder(para);
	for(int i=0;i<seqLen-1;i++){
		BRNode* nodeA = nodeList[i];
		BRNode* nodeB = nodeList[i+1];
		nodeA->connectToNeighbor = true;
		pb->buildPhosphate(nodeA->riboseConf, nodeB->riboseConf, nodeA->phoConf);
	}

	bool* connectToNeighbor = new bool[seqLen];
	for(int i=0;i<seqLen;i++){
		connectToNeighbor[i] = true;
	}
	connectToNeighbor[seqLen-1] = false;
	BRTreeInfo* info = new BRTreeInfo(seqLen, typeList, connectToNeighbor, nodeList, 0.0, rotLib);
	info->printPDB(output);

}




