/*
 * TestBasePairEnergy.cpp
 *
 *  Created on: 2023��8��23��
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
using namespace NSPpredNA;

int main(int argc, char** argv){
	clock_t start = clock();

	string pdbFile = string(argv[1]);
	string output = string(argv[2]);

	ofstream out;
	out.open(output.c_str(), ios::out);

	char xx[200];

	RNAPDB pdb = RNAPDB(pdbFile, "xxxx");


	
	vector<RNABase*> rawBaseList = pdb.getBaseList();

	vector<RNABase*> baseList;
	for(int i=0;i<rawBaseList.size();i++){
		RNABase* base = rawBaseList[i];
		if(base->baseTypeInt < 0 || base->baseTypeInt > 3) continue;
		base->updateCoordSystem();
		if(!base->hasLocalFrame) continue;
		baseList.push_back(rawBaseList[i]);
	}

	cout << "init et" << endl;
	ForceFieldPara* ffp = new ForceFieldPara();

	RnaEnergyTable* et = new RnaEnergyTable();
	et->loadEnergyWithout6D();

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
	BasePairLib* bpLib = new BasePairLib("xtb");

	for(int i=0;i<nodeList.size();i++){
		nodeA = nodeList[i];
		for(int j=i+1;j<nodeList.size();j++){
			sep = seqSep[j] - seqSep[i];
			if(sep > 2) sep = 2;
			nodeB = nodeList[j];

			BaseDistanceMatrix dm(nodeA->baseConf->cs1, nodeB->baseConf->cs1);
			int pairType = bpLib->getPairType(dm, nodeA->baseType, nodeB->baseType, 2);

			if(pairType < 0) continue;
			double eBB;
			
			eBB = bpLib->nnbEnergy[nodeA->baseType*4+nodeB->baseType][pairType];	

			//if(eBB > -1.6) continue;

			double eBBClash = baseBaseClash(nodeA, nodeB, sep, et, false);

			double eBR = getBaseRiboseEnergy(nodeA, nodeB, sep, et, false);
			eBR += getBaseRiboseEnergy(nodeB, nodeA, -sep, et, false);
			if(eBR > 0) eBR = 0;

			double eBP = getBasePhoEnergy(nodeA, nodeB, sep, et, false);
			eBP += getBasePhoEnergy(nodeB, nodeA, -sep, et, false);

			if(eBP > 0)
			 	eBP = 0;
			double eRR = getRiboseRiboseEnergy(nodeA, nodeB, sep, et, false);
			if(eRR > 0)
				eRR = 0;
			double eRP = getRibosePhoEnergy(nodeA, nodeB, sep, et, false);
			if(eRP > 0)
				eRP = 0;
			eRP += getRibosePhoEnergy(nodeB, nodeA, -sep, et, false);
			if(eRP > 0)
				eRP = 0;
			double ePP = getPhoPhoEnergy(nodeA, nodeB, sep, et, false);
			if(ePP > 0)
				ePP = 0;

			double eTot = eBB + eBR + eBP + eRR + eRP + ePP;
			if(eTot > 0)
				eTot = 0.0;
			
			
			sprintf(xx, "%d %d %d %4d %8.3f", nodeA->baseType, nodeB->baseType, sep, pairType, eTot);
			out << string(xx) << endl;
			
		}
	}

	out.close();


}

