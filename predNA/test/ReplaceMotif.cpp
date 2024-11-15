/*
 * ReplaceMotif.cpp
 *
 *  Created on: 2022��5��26��
 *      Author: pengx
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
		cout << "Usage: replaceMotif $INPUTFILE $OUTPUT" << endl;
		exit(0);
	}


	string inputFile = string(argv[1]);
	string output = string(argv[2]);

	RotamerLib* rotLib = new RotamerLib();

	NSPtools::InputParser input(inputFile);
	string pdbFile = input.getValue("pdb");
	string pdbFileM = input.getValue("pdbM");
	string baseSeq = input.getValue("seq");
	string fgA = input.getValue("fgA");
	string fgB = input.getValue("fgB");
	string mtf = input.getValue("mtf");

	int seqLen = baseSeq.length();
	BRNode** nodeList = new BRNode*[seqLen];

	RNAPDB pdb(pdbFile, "a");
	RNAPDB pdbM(pdbFileM, "b");
	vector<RNABase*> baseListP = pdb.getBaseList();
	vector<RNABase*> baseListM = pdbM.getBaseList();

	int motifLen = baseListM.size();

	cout << "parse fragment" << endl;

	if(fgA.length() != baseSeq.length()){
		cout << "fgA length error" << endl;
	}
	if(fgB.length() != baseSeq.length()){
		cout << "fgB length error" << endl;
	}

	int nA = 0;
	int nB = 0;

	int posX=-1;
	int posY=-1;
	map<int, int> seqPosToMtfPos;

	for(int i=0;i<baseSeq.length();i++){
		seqPosToMtfPos[i] = -1;
	}
	int m = 0;
	for(int i=0;i<baseSeq.length();i++){
		char c = mtf[i];
		if(c == '-')
			continue;

		seqPosToMtfPos[i] = m;
		m++;
	}
	if(m != motifLen) {
		cout << "motif length not match: " << motifLen << " " << m << endl;
	}

	for(int i=0;i<baseSeq.length();i++){
		if(fgA[i] == 'A')
			nA++;
		else if(fgA[i] != '-') {
			cout << "invalid fgA: " << fgA << endl;
			exit(0);
		}

		if(fgB[i] == 'B')
			nB++;
		else if(fgB[i] != '-') {
			cout << "invalid fgB: " << fgB << endl;
			exit(0);
		}

		if((fgA[i] == '-' && fgB[i] == '-') ||(fgA[i] == 'A' && fgB[i] == 'B') ){
			cout << "invalid fgA: " << fgA << endl;
			cout << "invalid fgB: " << fgB << endl;
			exit(0);
		}

		char c = mtf[i];
		if(c == 'X') {
			posX = i;
		}
		if(c == 'Y') {
			posY = i;
		}
	}

	if(posX < 0 || posY <0) {
		cout << "invalid mtf " << mtf << endl;
		exit(0);
	}


	map<char, int> c2i;
	c2i['A'] = 0;
	c2i['U'] = 1;
	c2i['G'] = 2;
	c2i['C'] = 3;
	int* typeList = new int[seqLen];
	LocalFrame csListP[seqLen];
	LocalFrame csListM[motifLen];

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
	for(int i=0;i<motifLen;i++){
		LocalFrame cs = baseListM[i]->getCoordSystem();
		csListM[i] = cs;
	}

	for(int i=0;i<seqLen;i++){
		if(fgA[i] == '-' && mtf[i] == '-'){
			LocalFrame csOld = csListP[i];
			CsMove mvXY = csListM[seqPosToMtfPos[posY]] - csListM[seqPosToMtfPos[posX]];
			CsMove mvYB = csListP[i] - csListP[posY];
			LocalFrame csNew = csListP[posX] + mvXY + mvYB;
			csListA[i] = csNew;
		}
		else if(fgA[i] == '-') {
			LocalFrame csOld = csListP[i];
			CsMove mvXM = csListM[seqPosToMtfPos[i]] - csListM[seqPosToMtfPos[posX]];
			LocalFrame csNew = csListP[posX] + mvXM;
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

