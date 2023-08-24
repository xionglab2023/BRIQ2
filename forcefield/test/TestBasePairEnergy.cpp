/*
 * TestBasePairEnergy.cpp
 *
 *  Created on: 2023Äê8ÔÂ23ÈÕ
 *      Author: nuc
 */



#include "geometry/localframe.h"
#include "geometry/TransMatrix.h"
#include "model/BaseDistanceMatrix.h"
#include "model/BasePairComposition.h"
#include "model/RiboseRotamerLib.h"
#include "model/AtomLib.h"
#include "predNA/BRNode.h"
#include <time.h>
#include <stdio.h>
#include "predNA/BRFoldingTree.h"
#include "model/StructureModel.h"
#include "forcefield/AtomicClashEnergy.h"
#include "forcefield/HbondEnergy.h"

#include "model/BasePairLib.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredna;

int main(int argc, char** argv){
	clock_t start = clock();

	string pdbFile = string(argv[1]);
	RNAPDB pdb = RNAPDB(pdbFile, "xxxx");
	RNAChain* rc = pdb.getFirstChain();
	vector<RNABase*> baseList = rc->getBaseList();
	cout << "start" << endl;
	XPara* para = new XPara();
	cout << "init et" << endl;
	ForceFieldPara* ffp = new ForceFieldPara();

	RnaEnergyTable* et = new RnaEnergyTable();
	AtomicClashEnergy* acET = new AtomicClashEnergy(ffp);
	cout << "finish init et" << endl;
	double tot = 0.0;
	double e;

	cout << "init rotLib" << endl;
	RotamerLib* rotLib = new RotamerLib();

	vector<BRNode*> nodeList;

	cout << "init seqSep" << endl;
	int seqSep[baseList.size()];
	for(int i=0;i<baseList.size();i++){
		seqSep[i] = 0;
	}

	cout << "init nodeList" << endl;

	int sep = 0;
	for(int i=0;i<baseList.size();i++){
		RNABase* baseA = baseList[i];
		seqSep[i] = sep;
		RiboseRotamer* riboRot;
		if(baseA->backboneComplete())
			riboRot = rotLib->riboseRotLib->getNearestRotamer(baseA);
		else
			riboRot = rotLib->riboseRotLib->getLowestEnergyRotamer(baseA->baseTypeInt);
		PhosphateRotamer* phoRot;
		if(i<baseList.size()-1 && baseList[i]->connectToNeighbor(baseList[i+1])){
			phoRot = new PhosphateRotamer(baseList[i], baseList[i+1]);
			sep ++;
		}
		else {
			phoRot = rotLib->phoRotLib->prLib[0][0];
			sep += 5;
		}
		BRNode* node = new BRNode(baseA, riboRot, phoRot, rotLib);
		nodeList.push_back(node);
	}

	BRNode* nodeA;
	BRNode* nodeB;
	int nA, nB;
	double dd;
	double clashEnergy;

	cout << "init basePair lib" << endl;
	BasePairLib* bpLib = new BasePairLib();

	cout << "calculate energy" << endl;
	for(int i=0;i<nodeList.size();i++){
		nodeA = nodeList[i];
		for(int j=i+1;j<nodeList.size();j++){
			sep = seqSep[j] - seqSep[i];
			if(sep > 3) sep = 3;
			nodeB = nodeList[j];


			double ePair = bpLib->getPairEnergy(baseList[i], baseList[j]);

			double eBB = getBaseBaseEnergy(nodeA, nodeB, sep, et, false);
			double eBBClash = baseBaseClash(nodeA, nodeB, sep, et, 0.0, false);
			double eBR = getBaseRiboseEnergy(nodeA, nodeB, sep, et, false);
			eBR += getBaseRiboseEnergy(nodeB, nodeA, -sep, et, false);
			double eBP = getBasePhoEnergy(nodeA, nodeB, sep, et, false);
			eBP += getBasePhoEnergy(nodeB, nodeA, -sep, et, false);
			double eRR = getRiboseRiboseEnergy(nodeA, nodeB, sep, et, false);
			double eRP = getRibosePhoEnergy(nodeA, nodeB, sep, et, false);
			eRP += getRibosePhoEnergy(nodeB, nodeA, -sep, et, false);
			double ePP = getPhoPhoEnergy(nodeA, nodeB, sep, et, false);

			double eTot = eBB + eBR + eBP + eRR + eRP + ePP;

			if(abs(ePair) > 0.1 || abs(eTot) > 0.1)
			{
				//printf("BB: %8.3f BR: %8.3f BP: %8.3f RR: %8.3f RP: %8.3f PP: %8.3f\n", eBB, eBR, eBP, eRR, eRP, ePP);
				printf("baseA: %3s baseB: %3s idA: %2d idB: %2d ePair: %8.3f eTot: %8.3f\n", baseList[i]->baseID.c_str(), baseList[j]->baseID.c_str(),i, j, ePair, eTot);
			}
		}
	}


}

