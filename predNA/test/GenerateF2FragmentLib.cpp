/*
 * GenerateF2FragmentLib.cpp
 *
 */

#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "forcefield/RnaAtomicEnergyTable.h"
#include "predNA/BaseMoveLibrary.h"
#include "predNA/BRNode.h"
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

	string moveType = string(argv[1]);

	RotamerLib* rotLib = new RotamerLib();
	XPara para;

	RnaAtomicEnergyTable* et = new RnaAtomicEnergyTable(&para);
	RiboseOxygenEnergyTable* rET = new RiboseOxygenEnergyTable();
	PO3Builder* pb = new PO3Builder(&para);
	ofstream out;
	ifstream in;
	string s;
	vector<string> spt;
	XYZ* tListB = new XYZ[11];
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			char xx[20];
			sprintf(xx, "%d", i*4+j);
			string inputfile = "/export/home/s2982206/aspen/rnaModeling/fragLib/"+moveType+"/"+moveType+".rot-"+string(xx);
			in.open(inputfile, ios::in);
			string outfile = "/export/home/s2982206/aspen/rnaModeling/fragLib/f2/"+moveType+"-"+string(xx)+".frag";
			out.open(outfile.c_str(), ios::out);
			getline(in, s);
			while(getline(in,s)){
				CsMove mv(s.substr(4, 131));
				LocalFrame csA;
				LocalFrame csB = csA + mv;
				BRNode* nodeA = new BRNode(i, 1, rotLib);
				nodeA->baseConf->updateCoords(csA);

				BRNode* nodeB = new BRNode(j, 2, rotLib);
				nodeB->baseConf->updateCoords(csB);

				splitString(s, " ", &spt);
				int idA = atoi(spt[14].c_str());
				int idB = atoi(spt[15].c_str());
				RiboseRotamer* rotA = rotLib->riboseRotLib->rotLib[i][idA];
				RiboseRotamer* rotB = rotLib->riboseRotLib->rotLib[j][idB];
				nodeA->riboseConf->updateRotamer(rotA);
				nodeB->riboseConf->updateRotamer(rotB);

				pb->buildPhosphate(nodeA->riboseConf, nodeB->riboseConf, nodeA->phoConf);

				char xxx[1000];
				sprintf(xxx, "%s %3d %3d %6.1f %6.1f\n", mv.toString().c_str(), idA, idB, nodeA->phoConf->rot->dihed1, nodeA->phoConf->rot->dihed2);
				out << string(xxx);

			}
			out.close();
			in.close();
		}
	}
}


