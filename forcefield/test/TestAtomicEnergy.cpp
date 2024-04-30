/*
 * TestAtomicEnergy.cpp
 *
 */


#include "geometry/localframe.h"
#include "geometry/TransMatrix.h"
#include "model/BaseDistanceMatrix.h"
#include "model/BasePairComposition.h"
#include "model/RiboseRotamerLib.h"
#include "forcefield/XPara.h"
#include "model/AtomLib.h"
#include "predNA/BRNode.h"
#include <time.h>
#include <stdio.h>

#include "model/StructureModel.h"
#include "forcefield/AtomicClashEnergy.h"
#include "forcefield/HbondEnergy.h"

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


	cout << "init atomLib" << endl;
	AtomLib* atLib = new AtomLib();

	cout << "calculate clash energy" << endl;
	for(int i=0;i<nodeList.size();i++){
		nodeA = nodeList[i];
		vector<string> nameA;
		atLib->getRnaSidechainAtoms(nodeA->baseType, nameA);
		for(int j=i+1;j<nodeList.size();j++){
			nodeB = nodeList[j];
			vector<string> nameB;
			atLib->getRnaSidechainAtoms(nodeB->baseType, nameB);
			int sep = seqSep[j] - seqSep[i];
			if(squareDistance(nodeA->baseConf->coords[0], nodeB->baseConf->coords[0]) < 256.0) {
				nA = nodeA->baseConf->rot->atomNum;
				nB = nodeB->baseConf->rot->atomNum;
				for(int k=0;k<nA;k++){
					for(int l=0;l<nB;l++){
						dd = squareDistance(nodeA->baseConf->coords[k], nodeB->baseConf->coords[l]);
						if(dd < 16.0) {
							clashEnergy = acET->getBaseBaseEnergy(nodeA->baseType, k, nodeB->baseType, l, dd, sep);
							if(abs(clashEnergy) > 0.1){
								printf("baseA: %3s %c baseB: %3s %c atomA: %3s atomB: %3s distance: %5.3f energy: %7.3f\n", baseList[i]->baseID.c_str(), baseList[i]->baseType, baseList[j]->baseID.c_str(), baseList[j]->baseType, nameA.at(k).c_str(), nameB.at(l).c_str(), sqrt(dd), clashEnergy);
							}
						}
					}
				}
			}
		}
	}


	HbondEnergy* hbET = new HbondEnergy(ffp);
	cout << "calculate hbond energy" << endl;

	for(int i=0;i<nodeList.size();i++){
		nodeA = nodeList[i];
		int polarNumA = nodeA->baseConf->rot->polarAtomNum;
		vector<string> nameA;
		atLib->getRnaSidechainAtoms(nodeA->baseType, nameA);
		for(int j=i+1;j<nodeList.size();j++){
			nodeB = nodeList[j];
			int polarNumB = nodeB->baseConf->rot->polarAtomNum;
			vector<string> nameB;
			atLib->getRnaSidechainAtoms(nodeB->baseType, nameB);
			if(squareDistance(nodeA->baseConf->coords[0], nodeB->baseConf->coords[0]) < 256.0) {
				for(int k=0;k<polarNumA;k++){
					int indexA = nodeA->baseConf->rot->polarAtomIndex[k];
					for(int l=0;l<polarNumB;l++){
						int indexB = nodeB->baseConf->rot->polarAtomIndex[l];
						double d = nodeA->baseConf->csPolar[k].origin_.distance(nodeB->baseConf->csPolar[l].origin_);
						double hbEne = hbET->getEnergy(nodeA->baseConf->rot->polarAtomUniqueID[k], nodeA->baseConf->csPolar[k], nodeB->baseConf->rot->polarAtomUniqueID[l], nodeB->baseConf->csPolar[l]);
						if(abs(hbEne)>0.2)
							printf("baseA: %3s %c baseB: %3s %c atomA: %3s atomB: %3s dist: %5.3f energy: %7.3f\n", baseList[i]->baseID.c_str(), baseList[i]->baseType, baseList[j]->baseID.c_str(), baseList[j]->baseType, nameA.at(indexA).c_str(), nameB.at(indexB).c_str(),d, hbEne);
					}
				}
			}
		}
	}



}
