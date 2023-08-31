/*
 * RNAScoring.cpp
 *
 */


#include "geometry/localframe.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

#include "forcefield/RnaEnergyTable.h"
#include "model/StructureModel.h"
#include "predNA/BRFoldingTree.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredna;
using namespace std;

int main(int argc, char** argv){
	string templatePDB = string(argv[1]);
	string pdbFile = string(argv[2]);
	string paraFile = string(argv[3]);

	RnaEnergyTable et(paraFile);
	//RiboConnectToPO3 et;

	RNAPDB pdb0(templatePDB, "xxxx");
	vector<RNABase*> baseList0 = pdb0.getBaseList();
	int seqLen = baseList0.size();
	bool connectToDownstream[seqLen];
	for(int i=0;i<seqLen;i++){
		connectToDownstream[i] = false;
	}
	for(int i=0;i<seqLen-1;i++){
		if(baseList0[i]->connectToNeighbor(baseList0[i+1]))
			connectToDownstream[i] = true;
	}
	RotamerLib* rotLib = new RotamerLib();



	{
		RNAPDB pdb(pdbFile, "xxxx");
		vector<RNABase*> baseList = pdb.getBaseList();

		cout << "PDB base number: " << baseList.size() << endl;

		if(baseList.size() != seqLen){
			cout << "RNA length not equal to template: " << pdbFile << endl;
		}
		vector<BRNode*> nodes;
		for(int i=0;i<seqLen;i++){
			RiboseRotamer* rot = rotLib->riboseRotLib->getNearestRotamer(baseList[i]);
			BRNode* node = new BRNode(baseList[i], rotLib);
			node->riboseConf->updateRotamer(rot);
			node->riboseConfTmp->updateRotamer(rot);
			node->connectToNeighbor = connectToDownstream[i];
			nodes.push_back(node);
		}

		for(int i=0;i<seqLen-1;i++){
			BRNode* nodeA = nodes[i];
			BRNode* nodeB = nodes[i+1];
			if(!connectToDownstream[i]) continue;
			et.pb->buildPhosphate(nodeA->riboseConf, nodeB->riboseConf, nodeA->phoConf);
		}


		double bpEnergy1=0;
		double bpEnergy2=0;
		double bpEnergy3=0;

		double baseClash = 0;
		double oxyE = 0;

		double phoConnectionEnergy = 0;
		double rotEnergy = 0;

		for(int i=0;i<seqLen;i++){
			rotEnergy += nodes[i]->riboseConf->rot->energy;
			phoConnectionEnergy += nodes[i]->phoConf->ene;


		}

		for(int i=0;i<seqLen;i++){
			cout << "node: " << i << " rot: " << nodes[i]->riboseConf->rot->energy << endl;
		}

		for(int i=0;i<seqLen;i++){
			cout << "node: " << i << " pho: " << nodes[i]->phoConf->ene << endl;
		}

		double shift = 0;

		for(int a=0;a<seqLen;a++){
			BRNode* nodeA = nodes[a];

			for(int b=a+1;b<seqLen;b++){
				BRNode* nodeB = nodes[b];
				int sep = 3;
				if(b == a+1 && connectToDownstream[a])
					sep = 1;
				if(b == a+2 && connectToDownstream[a] && connectToDownstream[a+1])
					sep = 2;

				int sep2 = -sep;
				if(sep2 == -3) sep2 = 3;

				int i,j;

				double clash = 0;
				double bb = 0;
				double eo = 0;

				double minDD = 990009.9;


				/*
				 * base base energy
				 */
				baseClash += baseBaseClash(nodeA, nodeB, sep, &et,false);
				clash = baseBaseClash(nodeA, nodeB, sep, &et, false);

				double e1,e2,e3,e4,e5,e6,e7,e8;
				e1 = getBaseRiboseEnergy(nodeA, nodeB, sep, &et, false);
				e2 = getBaseRiboseEnergy(nodeB, nodeA, sep2, &et, false);
				e3 = getBasePhoEnergy(nodeA, nodeB, sep, &et,false);
				e4 = getBasePhoEnergy(nodeB, nodeA, sep2, &et, false);
				e5 = getRiboseRiboseEnergy(nodeA, nodeB, sep, &et, false);
				e6 = getRibosePhoEnergy(nodeA, nodeB, sep, &et, false);
				e7 = getRibosePhoEnergy(nodeB, nodeA, sep2, &et, false);
				e8 = getPhoPhoEnergy(nodeA, nodeB, sep, &et, false);

				if(abs(e1) > 0.05){
					printf("e1: %7.3f\n", e1);
				}
				if(abs(e2) > 0.05){
					printf("e1: %7.3f\n", e2);
				}
				if(abs(e3) > 0.05){
					printf("e1: %7.3f\n", e3);
				}
				if(abs(e4) > 0.05){
					printf("e1: %7.3f\n", e4);
				}
				if(abs(e5) > 0.05){
					printf("e1: %7.3f\n", e5);
				}
				if(abs(e6) > 0.05){
					printf("e1: %7.3f\n", e6);
				}
				if(abs(e7) > 0.05){
					printf("e1: %7.3f\n", e7);
				}
				if(abs(e8) > 0.05){
					printf("e1: %7.3f\n", e8);
				}


				oxyE += eo;


				if(sep == 1){
					//bp energy
					bb = getBaseBaseEnergy(nodeA, nodeB, sep, &et, false);
					bpEnergy1 += bb;
				}
				else if(sep == 2){
					//bp energy
					bb = getBaseBaseEnergy(nodeA, nodeB, sep, &et, false);
					bpEnergy2 += bb;
					}
				else if(sep == 3){
					//bp energy
					bb = getBaseBaseEnergy(nodeA, nodeB, sep, &et, false);
					bpEnergy3 += bb;
				}

				double t = bb + eo + clash;

				if(abs(t) > 0.05)
					printf("node: %2d node: %2d baseBaseE: %7.3f oxyE: %7.3f clash: %7.3f\n", a, b, bb, eo, clash);
			}
		}

		//double tot = bpEnergy1 + bpEnergy2 + bpEnergy3 + oxyE + phoConnectionEnergy + rotEnergy + baseClash;
		//printf("E: %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n", bpEnergy1, bpEnergy2, bpEnergy3, oxyE, phoConnectionEnergy, rotEnergy, baseClash, tot);


		for(int i=0;i<seqLen;i++){
			delete nodes[i];
		}
	}



}



