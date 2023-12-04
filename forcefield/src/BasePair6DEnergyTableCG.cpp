/*
 * BasePair6DEnergyTableCG.cpp
 *
 *  Created on: 2023Äê11ÔÂ28ÈÕ
 *      Author: nuc
 */

#include "forcefield/BasePair6DEnergyTableCG.h"

namespace NSPforcefield {

BasePair6DEnergyTableCG::BasePair6DEnergyTableCG(ForceFieldPara* para) {

	this->wtNb = para->wtBp1;
	this->wtNnb = para->wtBp2;

	string path = NSPdataio::datapath();

	int indexA, indexB, clusterID;
	double ene;
	ifstream file;
	string augc = "AUGC";

	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			string pairType = augc.substr(i,1) + augc.substr(j,1);
			file.open(path + "pairEneCG/nb/"+pairType+".ene");

			if(!file.is_open()) {
				cout << "can't open file " << path + "pairEne/nb/" +pairType+".ene" << endl;
			}
			while(file >> indexA >> indexB >> ene >> clusterID){
				this->nbKeysEnergy[(i*4+j)*2250+indexA][indexB] = ene;
			}
			file.close();

			file.open(path + "pairEneCG/nnb/"+pairType+".ene");
			if(!file.is_open()) {
				cout << "can't open file " << path + "pairEne/nnb/" +pairType+".ene" << endl;
			}
			while(file >> indexA >> indexB >> ene >> clusterID){
				this->nnbKeysEnergy[(i*4+j)*2250+indexA][indexB] = ene;
			}
			file.close();

		}
	}
}

BasePair6DEnergyTableCG::~BasePair6DEnergyTableCG() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPpredNA */
