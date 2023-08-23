/*
 * GenerateF3FragmentLib.cpp
 *
 */

#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "forcefield/RnaAtomicEnergyTable.h"
#include "predNA/BaseMoveLibrary.h"
#include "predNA/BRNode.h"
#include "predNA/ThreeBaseMoveLibrary.h"
#include "model/BaseRotamerLib.h"
#include "model/RiboseRotamerLib.h"
#include "forcefield/PO3Builder.h"
#include "forcefield/RiboseOxygenEnergyTable.h"


using namespace NSPpredna;
using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPgeometry;

double getBaseRiboseEnergy(BRNode* nodeA, BRNode* nodeB, RnaAtomicEnergyTable* et, int sep){
	double clashEnergy = 0.0;
	int i,j;
	for(i=0;i<nodeA->baseConf->rot->atomNum;i++){
		for(j=0;j<8;j++){
			clashEnergy += et->getBaseRiboseEnergy(nodeA->baseType, i, j, squareDistance(nodeA->baseConf->coords[i], nodeB->riboseConf->coords[j]), sep);
		}
	}
	return clashEnergy;
}

double getBasePhoEnergy(BRNode* nodeA, BRNode* nodeB, RnaAtomicEnergyTable* et, int sep){
	double clashEnergy = 0.0;
	int i,j;
	for(i=0;i<nodeA->baseConf->rot->atomNum;i++){
		for(j=1;j<4;j++){
			clashEnergy += et->getBasePhoEnergy(nodeA->baseType, i, j, squareDistance(nodeA->baseConf->coords[i], nodeB->phoConf->coords[j]), sep);
		}
	}
	return clashEnergy;
}

double getRiboseRiboseEnergy(BRNode* nodeA, BRNode* nodeB, RnaAtomicEnergyTable* et, int sep){
	double clashEnergy = 0.0;
	int i,j;
	for(i=0;i<8;i++){
		for(j=0;j<8;j++){
			clashEnergy += et->getRiboseRiboseEnergy(i,j,squareDistance(nodeA->riboseConf->coords[i], nodeB->riboseConf->coords[j]), sep);
		}
	}
	return clashEnergy;
}

double getRibosePhoEnergy(BRNode* nodeA, BRNode* nodeB, RnaAtomicEnergyTable* et, int sep){
	int i,j;
	double clashEnergy = 0;
	for(i=0;i<8;i++){
		for(j=1;j<4;j++){
			clashEnergy += et->getRibosePhoEnergy(i,j,squareDistance(nodeA->riboseConf->coords[i], nodeB->phoConf->coords[j]), sep);
		}
	}
	return clashEnergy;
}

double getPhoPhoEnergy(BRNode* nodeA, BRNode* nodeB, RnaAtomicEnergyTable* et, int sep){
	int i,j;
	double clashEnergy = 0;
	for(i=1;i<4;i++){
		for(j=1;j<4;j++){
			clashEnergy += et->getPhoPhoEnergy(i,j,squareDistance(nodeA->phoConf->coords[i], nodeB->phoConf->coords[j]), sep);
		}
	}
	return clashEnergy;
}



int main(int argc, char** argv){

	string outpath = string(argv[1]);


	RotamerLib* rotLib = new RotamerLib();
	XPara* para = new XPara();
	RnaAtomicEnergyTable* et = new RnaAtomicEnergyTable(para);
	RiboseOxygenEnergyTable* rET = new RiboseOxygenEnergyTable();
	PO3Builder* pb = new PO3Builder(para);

	ofstream out;
	ifstream in;
	string s;
	vector<string> spt;

	XYZ* tListB = new XYZ[11];
	XYZ* tListC = new XYZ[11];
	string augc = "AUGC";


	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			for(int k=0;k<4;k++){
				string typeTag = ""+augc.substr(i,1)+augc.substr(j,1)+augc.substr(k,1);

				char xx[20];
				sprintf(xx, "%d", i*16+j*4+k);
				string outfile = outpath + "/f3-all-" + string(xx) + ".frag";

				string inputfile = "/export/home/s2982206/aspen/rnaModeling/fragLib/f3All/"+typeTag+".rot";
				in.open(inputfile, ios::in);
				out.open(outfile.c_str(), ios::out);

				getline(in, s);
				while(getline(in,s)){
					CsMove mv1(s.substr(6, 137)); //cm12
					CsMove mv2(s.substr(138, 268)); //cm23
					CsMove mv3(s.substr(270, 400)); //cm13
					LocalFrame csA;
					LocalFrame csB = csA + mv1;
					LocalFrame csC = csA + mv3;
					BRNode* nodeA = new BRNode(i, 1, rotLib);
					nodeA->baseConf->updateCoords(csA);

					BRNode* nodeB = new BRNode(j, 2, rotLib);
					nodeB->baseConf->updateCoords(csB);

					BRNode* nodeC = new BRNode(j, 2, rotLib);
					nodeC->baseConf->updateCoords(csC);

					splitString(s, " ", &spt);
					int idA = atoi(spt[39].c_str());
					int idB = atoi(spt[40].c_str());
					int idC = atoi(spt[41].c_str());

					RiboseRotamer* rotA = rotLib->riboseRotLib->rotLib[i][idA];
					RiboseRotamer* rotB = rotLib->riboseRotLib->rotLib[j][idB];
					RiboseRotamer* rotC = rotLib->riboseRotLib->rotLib[k][idC];

					nodeA->riboseConf->updateRotamer(rotA);
					nodeB->riboseConf->updateRotamer(rotB);
					nodeC->riboseConf->updateRotamer(rotC);

					nodeA->riboseConf->updateLocalFrame(csA);
					nodeB->riboseConf->updateLocalFrame(csB);
					nodeC->riboseConf->updateLocalFrame(csC);

					pb->buildPhosphate(nodeA->riboseConf, nodeB->riboseConf, nodeA->phoConf);
					pb->buildPhosphate(nodeB->riboseConf, nodeC->riboseConf, nodeB->phoConf);

					char xxx[1000];
					//mv12, mv23, rotA, rotB, rotC, phoA, phoB
					sprintf(xxx, "%s %s %3d %3d %3d %6.1f %6.1f %6.1f %6.1f\n", mv1.toString().c_str(), mv2.toString().c_str(), idA, idB, idC, nodeA->phoConf->rot->dihed1, nodeA->phoConf->rot->dihed2, nodeB->phoConf->rot->dihed1, nodeB->phoConf->rot->dihed2);
					out << string(xxx);
				}
				in.close();
				out.close();
			}
		}
	}
}


