/*
 * AddBaseToTarget.cpp
 *
 *  Created on: 2022��5��25��
 *      Author: pengx
 */


#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include "model/StructureModel.h"
#include "predNA/BRFoldingTree.h"
#include "predNA/MCRun.h"
#include "forcefield/ForceFieldPara.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredna;
using namespace std;

int main(int argc, char** argv){

	if(argc !=  3 || argv[1] == "-h"){
		cout << "Usage: add_base $INPUTFILE $OUTPUT" << endl;
		exit(0);
	}

	string inputFile = string(argv[1]);
	string output = string(argv[2]);

	RotamerLib* rotLib = new RotamerLib();

	NSPtools::InputParser input(inputFile);
	string pdbFileA = input.getValue("pdbA");
	string baseSeq = input.getValue("seq");
	vector<string> addList = input.getMultiValues("add");

	string fgA = input.getValue("fgA");


	cout << baseSeq << endl;
	cout << fgA << endl;

	int seqLen = baseSeq.length();
	BRNode** nodeList = new BRNode*[seqLen];

	RNAPDB pdbA(pdbFileA, "a");
	vector<RNABase*> baseListA = pdbA.getBaseList();

	cout << "parse fragment" << endl;

	if(fgA.length() != baseSeq.length()){
		cout << "fgA length error" << endl;
	}

	int nA = 0;
	int nB = 0;
	int nX = 0;
	int overlapPos = -1;
	int fragmentBNum = 0;
	for(int i=0;i<baseSeq.length();i++){
		if(fgA[i] == 'A')
			nA++;
		else if(fgA[i] != '-') {
			cout << "invalid fgA: " << fgA << endl;
			exit(0);
		}
	}

	if(nA != baseListA.size()) {
		cout << "fgA not equal to pdbA" << endl;
		exit(0);
	}

	map<char, int> c2i;
	c2i['A'] = 0;
	c2i['U'] = 1;
	c2i['G'] = 2;
	c2i['C'] = 3;
	int* typeList = new int[seqLen];
	LocalFrame csListA[seqLen];
	bool hasCs[seqLen];

	for(int i=0;i<fgA.length();i++){
		if(fgA[i] == 'A')
			hasCs[i] = true;
		else
			hasCs[i] = false;
	}

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
		char c = fgA[i];
		if(c == 'A'){
			baseListA[id]->updateCoordSystem();
			if(!baseListA[id]->hasLocalFrame){
				cout << "pdbA not complete: " << i  << endl;
				exit(0);
			}

			LocalFrame cs = baseListA[id]->getCoordSystem();
			baseList[i] = baseListA[id];
			id++;

			csListA[i] = cs;
		}
	}

	CsMove cmWcNb("1.179875   3.813016   -3.749098  0.8472664  0.5174403  -0.1199800 -0.5178250 0.8549461  0.0304035  0.1183084  0.0363688  0.9923107");
	CsMove cmWcNbRev = cmWcNb.reverse();
	CsMove cmWc("5.629524   8.758814   0.586972   0.4021382  -0.9013801 -0.1606198 -0.8781193 -0.4293707 0.2110625  -0.2592130 0.0561671  -0.9641856");
	CsMove cmCrossNb("-6.045276  2.905033   -1.371541  -0.9439384 -0.3253931 -0.0556740 -0.2306567 0.5294329  0.8163935  -0.2361731 0.7834668  -0.5748061");
	CsMove cmCrossNbRev("-5.360224  -2.430554  -3.496585  -0.9439384 -0.2306567 -0.2361731 -0.3253931 0.5294329  0.7834668  -0.0556740 0.8163935  -0.5748061");

	vector<string> spt;
	for(int i=0;i<addList.size();i++){
		splitString(addList[i], " ", &spt);
		int idA = atoi(spt[0].c_str());
		int idB = atoi(spt[1].c_str());

		string ctType = spt[2];
		CsMove cm;
		if(ctType == "wc")
			cm = cmWc;
		else if(ctType == "wcNb")
			cm = cmWcNb;
		else if(ctType == "revWcNb")
			cm = cmWcNbRev;
		else if(ctType == "crossNb")
			cm = cmCrossNb;
		else if(ctType == "revCrossNb")
			cm = cmCrossNbRev;
		else {
			cout << "ct not support: " << ctType << endl;
			cout << "wc" << endl;
			cout << "wcNb" << endl;
			cout << "revWcNb" << endl;
		}


		if(!hasCs[idA]) {
			cout << "invalid action " << addList[i] << " base A not built" << endl;
		}
		if(hasCs[idB]) {
			cout << "invalid action " << addList[i] << " base B already built" << endl;
		}

		csListA[idB] = csListA[idA] + cm;
		hasCs[idB] = true;
	}

	for(int i=0;i<seqLen;i++){
		if(!hasCs[i]) {
			cout << "pos " << i << "not built" << endl;
		}
	}


	cout << "generate nodes" << endl;
	for(int i=0;i<seqLen;i++){
		cout << "node: " << i << endl;
		BRNode* node = new BRNode(typeList[i], i, rotLib);
		node->baseConf->updateCoords(csListA[i]);

		RiboseRotamer* rot;

		if(fgA[i] == 'A') {
			if(!baseList[i]->backboneComplete())
				rot = rotLib->riboseRotLib->getLowestEnergyRotamer(baseList[i]->baseTypeInt);
			else
				rot = new RiboseRotamer(baseList[i]);
		}
		else {
			rot = rotLib->riboseRotLib->getLowestEnergyRotamer(typeList[i]);
		}

		node->riboseConf->updateLocalFrame(csListA[i]);
		node->riboseConf->updateRotamer(rot);

		nodeList[i] = node;
	}

	cout << "build po3" << endl;

	ForceFieldPara* para = new ForceFieldPara();

	RiboseOxygenEnergyTable* rET = new RiboseOxygenEnergyTable(para);
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



