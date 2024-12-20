/*
 * AssembleTwoFragments.cpp
 *
 *  Created on: 2022��5��12��
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

	if(argc !=  3 || argv[1] == "-h"){
		cout << "Usage: rna_assemble $INPUTFILE $OUTPUT" << endl;
		exit(0);
	}

	string inputFile = string(argv[1]);
	string output = string(argv[2]);

	cout << "init rotLib" << endl;
	RotamerLib* rotLib = new RotamerLib();
	AtomLib* atLib = new AtomLib();

	cout << "parse input" << endl;

	NSPtools::InputParser input(inputFile);
	string pdbFileA = input.getValue("pdbA");
	string pdbFileB = input.getValue("pdbB");
	string baseSeq = input.getValue("seq");
	string fgA = input.getValue("fgA");
	string fgB = input.getValue("fgB");

	int seqLen = baseSeq.length();

	NuNode** nodeList = new NuNode*[seqLen];

	cout << "init pdb" << endl;
	RNAPDB pdbA(pdbFileA, "a");
	RNAPDB pdbB(pdbFileB, "b");
	vector<RNABase*> baseListA = pdbA.getBaseList();
	vector<RNABase*> baseListB = pdbB.getBaseList();

	cout << "parse fragment" << endl;
	cout << fgA << endl;
	cout << fgB << endl;

	if(fgA.length() != baseSeq.length()){
		cout << "fgA length error" << endl;
	}
	if(fgB.length() != baseSeq.length()){
		cout << "fgB length error" << endl;
	}

	int nA = 0;
	int nB = 0;
	int nX = 0;
	int overlapPos = 0;
	int fragmentBNum = 0;
	for(int i=0;i<baseSeq.length();i++){
		if(fgA[i] == 'A')
			nA++;
		else if(fgA[i] != '-') {
			cout << "invalid fgA: " << fgA << endl;
			exit(0);
		}

		if(fgB[i] == 'B')
			nB++;
		else if(fgB[i] == 'X') {
			nB++;
			nX++;
			overlapPos = i;
		}
		else if(fgB[i] != '-') {
			cout << "invalid fgB: " << fgB << endl;
			exit(0);
		}

		if(fgA[i] == '-' && fgB[i] == '-'){
			cout << "invalid fgA: " << fgA << endl;
			cout << "invalid fgB: " << fgB << endl;
			exit(0);
		}
	}

	if(nX != 1 && nX != 0){
		cout << "invalid fgB: " << fgB << endl;
		exit(0);
	}

	if(nA != baseListA.size()) {
		cout << "fgA not equal to pdbA" << endl;
		exit(0);
	}

	if(nB != baseListB.size()) {
		cout << "fgB not equal to pdbB" << endl;
		exit(0);
	}

	map<char, int> c2i;
	c2i['A'] = 0;
	c2i['U'] = 1;
	c2i['G'] = 2;
	c2i['C'] = 3;
	int* typeList = new int[seqLen];
	LocalFrame csListA[seqLen];
	LocalFrame csListB[seqLen];
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

	cout << "read coordinate system" << endl;

	id = 0;
	for(int i=0;i<seqLen;i++){
		char c = fgB[i];
		if(c == 'B' || c == 'X'){
			baseListB[id]->updateCoordSystem();
			if(!baseListB[id]->hasLocalFrame){
				cout << "pdbB not complete: " << i  << endl;
				exit(0);
			}
			LocalFrame cs = baseListB[id]->getCoordSystem();
			if(fgA[i] == '-'){
				baseList[i] = baseListB[id];
			}
			id++;
			csListB[i] = cs;
		}
	}

	for(int i=0;i<seqLen;i++){
		char c = fgB[i];
		if(c == 'B' || c == 'X'){
			LocalFrame csB = csListB[i];
			CsMove cm = csB - csListB[overlapPos];
			LocalFrame csB2 = csListA[overlapPos] + cm;
			cout << "pos: " << i << endl;
			csListA[i] = csB2;
		}
	}

	cout << "generate nodes" << endl;
	for(int i=0;i<seqLen;i++){
		cout << "node: " << i << endl;
		
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
		BaseRotamer* baseRot = rotLib->baseRotLib->baseLib[typeList[i]];
		NuNode* node = new NuNode(i, typeList[i], csListA[i], baseRot, rot, atLib);
		nodeList[i] = node;
	}

	cout << "build po3" << endl;

	ForceFieldPara* para = new ForceFieldPara();
	PO3Builder* pb = new PO3Builder(para);
	for(int i=0;i<seqLen-1;i++){
		NuNode* nodeA = nodeList[i];
		NuNode* nodeB = nodeList[i+1];
		pb->buildPhosphate(nodeA->riboseConf, nodeB->riboseConf, nodeA->phoConf);
	}


	bool* connectToNeighbor = new bool[seqLen];
	bool* fixed = new bool[seqLen];
	for(int i=0;i<seqLen;i++){
		connectToNeighbor[i] = true;
		fixed[i] = false;
	}
	connectToNeighbor[seqLen-1] = false;

	graphInfo* info = new graphInfo(seqLen, typeList, connectToNeighbor, fixed, nodeList, 0.0, atLib, 0);
	info->printPDB(output);

}


