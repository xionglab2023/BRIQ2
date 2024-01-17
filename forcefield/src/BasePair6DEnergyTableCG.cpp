/*
 * BasePair6DEnergyTableCG.cpp
 *
 *  Created on: 2023��11��28��
 *      Author: nuc
 */

#include "forcefield/BasePair6DEnergyTableCG.h"

namespace NSPforcefield {

BasePair6DEnergyTableCG::BasePair6DEnergyTableCG(ForceFieldPara* para, bool withBinary, int binaryMode):
cm2Key(withBinary, binaryMode)
{

	this->wtNb = para->wtBp1;
	this->wtNnb = para->wtBp2;

	string path = NSPdataio::datapath();

	int indexA, indexB, clusterID;
	double ene;
	ifstream file;
	string augc = "AUGC";

	if(withBinary && binaryMode==2) {
		string fileName = path + "../binaryData/pairEneCG/nb";
		ifstream ins;
        ins.open(fileName,ios::in | ios::binary);
        if(!ins.is_open()) {
            throw("Unable to open " + fileName);
        }
        auto* bb = new BinaryBook;
        bb->read(ins);
        ins.close();
		for(int i=0; i<4; i++) {
			for(int j=0;j<4;j++) {
				string tN = string(augc.substr(i,1)) + augc.substr(j,1) + ".ene";
				BinaryTable* btab = bb->tables_map.at(tN);
				auto& indexACol = get<BinaryColumn<int>>(*(btab->cols[0]));
				auto& indexBCol = get<BinaryColumn<int>>(*(btab->cols[1]));
				auto& eneCol = get<BinaryColumn<double>>(*(btab->cols[2]));
				// nbKeysEnergy 在元素map中随机赋值，且各map size 未知，无法并行
				for(int k=0;k<btab->nRow;k++) {
					this->nbKeysEnergy[(i*4+j)*2250+indexACol[k]][indexBCol[k]] = eneCol[k];
				}
			}
		}
		delete bb;

		fileName = path + "../binaryData/pairEneCG/nnb";
        ins.open(fileName,ios::in | ios::binary);
        if(!ins.is_open()) {
            throw("Unable to open " + fileName);
        }
        bb = new BinaryBook;
        bb->read(ins);
        ins.close();
		for(int i=0; i<4; i++) {
			for(int j=0;j<4;j++) {
				string tN = string(augc.substr(i,1)) + augc.substr(j,1) + ".ene";
				BinaryTable* btab = bb->tables_map.at(tN);
				auto& indexACol = get<BinaryColumn<int>>(*(btab->cols[0]));
				auto& indexBCol = get<BinaryColumn<int>>(*(btab->cols[1]));
				auto& eneCol = get<BinaryColumn<double>>(*(btab->cols[2]));
				for(int k=0;k<btab->nRow;k++) {
					this->nnbKeysEnergy[(i*4+j)*2250+indexACol[k]][indexBCol[k]] = eneCol[k];
				}
			}
		}
		delete bb;		
	} else if(withBinary && binaryMode == 1) {

	} else {
		for(int i=0;i<4;i++){
			for(int j=0;j<4;j++){
				string pairType = augc.substr(i,1) + augc.substr(j,1);
				file.open(path + "pairEneCG/nb/"+pairType+".ene");

				if(!file.is_open()) {
					cout << "can't open file " << path + "pairEneCG/nb/" +pairType+".ene" << endl;
				}
				while(file >> indexA >> indexB >> ene >> clusterID){
					this->nbKeysEnergy[(i*4+j)*2250+indexA][indexB] = ene;
				}
				file.close();

				file.open(path + "pairEneCG/nnb/"+pairType+".ene");
				if(!file.is_open()) {
					cout << "can't open file " << path + "pairEneCG/nnb/" +pairType+".ene" << endl;
				}
				while(file >> indexA >> indexB >> ene >> clusterID){
					this->nnbKeysEnergy[(i*4+j)*2250+indexA][indexB] = ene;
				}
				file.close();

			}
		}
	}
}

BasePair6DEnergyTableCG::~BasePair6DEnergyTableCG() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPpredNA */
