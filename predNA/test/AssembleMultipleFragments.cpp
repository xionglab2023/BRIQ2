/*
 * AssembleMultipleFragments.cpp
 *
 *  Created on: 2022��6��14��
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
		cout << "Usage: rna_assemble_multiFragments $INPUTFILE $OUTPUT" << endl;
		exit(0);
	}

	string inputFile = string(argv[1]);
	string output = string(argv[2]);

	cout << "assemble multiple fragments" << endl;
	RotamerLib* rotLib = new RotamerLib();

	NSPtools::InputParser input(inputFile);
	int fragmentNum = atoi(input.getValue("N").c_str());
	string pdbFileA = input.getValue("pdbA");
	vector<string> fileList;
	char xx[20];
	for(int i=1;i<fragmentNum+1;i++){
		sprintf(xx, "pdb%d", i);
		string tag = string(xx);
		fileList.push_back(input.getValue(tag));
	}

	cout << "fragment num: " << fragmentNum << endl;

	cout << "get seq" << endl;
	string baseSeq = input.getValue("seq");
	cout << "get fgA" << endl;
	string fgA = input.getValue("fgA");
	cout << "get fgB" << endl;
	string fgB = input.getValue("fgB");

	int seqLen = baseSeq.length();
	BRNode** nodeList = new BRNode*[seqLen];


	cout << "fileList " << fileList.size() << endl;
	RNAPDB* pdbA = new RNAPDB(pdbFileA, "a");
	vector<RNAPDB*> pdbBList;
	for(int i=0;i<fragmentNum;i++){
		RNAPDB* pdb = new RNAPDB(fileList[i], "b");
		pdbBList.push_back(pdb);
	}

	cout << "baseList A" << endl;

	vector<RNABase*> baseListA = pdbA->getBaseList();

	cout << "baseList B" << endl;
	vector<vector<RNABase*>> baseListListB;
	for(int i=0;i<fragmentNum;i++){
		baseListListB.push_back(pdbBList[i]->getBaseList());
	}

	cout << "parse fragment" << endl;

	if(fgA.length() != baseSeq.length()){
		cout << "fgA length error" << endl;
	}
	if(fgB.length() != baseSeq.length()){
		cout << "fgB length error" << endl;
	}

	int nA = 0;
	vector<int> nBList;
	vector<int> nXList;
	vector<int> overlapList;
	for(int i=0;i<fragmentNum;i++){
		nBList.push_back(0);
		nXList.push_back(0);
		overlapList.push_back(0);
	}

	for(int i=0;i<baseSeq.length();i++){
		if(fgA[i] == 'A')
			nA++;
		else if(fgA[i] == 'X'){
			nA++;
			if(fgB[i] >= '1' && fgB[i] <= '0'+fragmentNum){
				nXList[fgB[i]-'1']++;
				overlapList[fgB[i]-'1'] = i;
			}
			else {
				cout << "X match to gap" << endl;
				cout << "invalid fgA: " << fgA << endl;
				cout << "invalid fgB: " << fgB << endl;
				exit(0);
			}
		}
		else if(fgA[i] != '-') {
			cout << "invalid fgA: " << fgA << endl;
			exit(0);
		}

		if(fgB[i] >= '1' && fgB[i] <= '0'+fragmentNum)
			nBList[fgB[i]-'1']++;
		else if(fgB[i] != '-') {
			cout << fgB[i] << endl;
			cout << "invalid fgB: " << fgB << endl;
			exit(0);
		}

		if(fgA[i] == '-' && fgB[i] == '-'){
			cout << "invalid fgA: " << fgA << endl;
			cout << "invalid fgB: " << fgB << endl;
			exit(0);
		}
	}

	for(int i=0;i<fragmentNum;i++){
		int nX = nXList[i];
		if(nX != 1 && nX != 0){
			cout << "invalid fgB: " << fgB << endl;
			exit(0);
		}
	}


	if(nA != baseListA.size()) {
		cout << "fgA not equal to pdbA" << endl;
		exit(0);
	}

	for(int i=0;i<fragmentNum;i++){
		int nB = nBList[i];
		if(nB != baseListListB[i].size()) {
			cout << "fg" << (i+1) << " not equal to pdb" << (i+1) << endl;
			exit(0);
		}
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
		if(c == 'A' || c == 'X'){
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

	for(int k=1;k<fragmentNum+1;k++){
		id = 0;
		char fragID = '0'+k;
		for(int i=0;i<seqLen;i++){
			char c = fgB[i];
			if(c == fragID){
				baseListListB[k-1][id]->updateCoordSystem();
				if(!baseListListB[k-1][id]->hasLocalFrame){
					cout << "pdbB not complete: " << i  << endl;
					exit(0);
				}
				LocalFrame cs = baseListListB[k-1][id]->getCoordSystem();
				if(fgA[i] == '-'){
					baseList[i] = baseListListB[k-1][id];
				}
				id++;
				csListB[i] = cs;
			}
		}
	}


	for(int i=0;i<seqLen;i++){
		char c = fgA[i];
		if(c == '-'){
			int fragID = fgB[i]-'1';
			LocalFrame csB = csListB[i];
			CsMove cm = csB - csListB[overlapList[fragID]];
			LocalFrame csB2 = csListA[overlapList[fragID]] + cm;
			csListA[i] = csB2;
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





