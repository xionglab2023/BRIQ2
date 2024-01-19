/*
 * BasePair6DEnergyTableCG.cpp
 *
 *  Created on: 2023��11��28��
 *      Author: nuc
 */

#include "forcefield/BasePair6DEnergyTableCG.h"
#include "dataio/dirOperations.h"
#include <string.h>

namespace NSPforcefield {

using namespace NSPdataio;

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
		return;
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

int BasePair6DEnergyTableCG::dump() {
	// serialized dump method
	string outpath = datapath() + "../binaryCache";
	string fileName = "BasePair6DEnergyTableCG";
	char z0[4] = {'\0'};
	if(! makeDirs(outpath)) {
		throw("[Error]Unable to create " + outpath);
	}
	ofstream outs;
	outs.open(outpath + "/" + fileName, ios::out|ios::binary);
	if(!outs.is_open()) {
		throw("[Error]Fail to open " + outpath + "/" + fileName);
	}
	cm2Key.dump(outs);
	outs.write(reinterpret_cast<char*>(&wtNb), sizeof(double));
	outs.write(reinterpret_cast<char*>(&wtNnb), sizeof(double));
	int nbMapSizes[36000], nnbMapSizes[36000];
	for(int i=0; i<36000; i++) {
		nbMapSizes[i] = nbKeysEnergy[i].size();
		nnbMapSizes[i] = nnbKeysEnergy[i].size();
	}
	outs.write(reinterpret_cast<char*>(nbMapSizes),sizeof(int)*36000);
	for(int i=0; i<36000; i++) {
		int keys[nbMapSizes[i]];
		double vals[nbMapSizes[i]];
		int j = 0;
		for(auto& it1:nbKeysEnergy[i]) {
			keys[j] = it1.first;
			vals[j] = it1.second;
			j++;
		}
		outs.write(reinterpret_cast<char*>(keys), sizeof(int)*nbMapSizes[i]);
		if(nbMapSizes[i]%2) {
			outs.write(reinterpret_cast<char*>(z0), sizeof(int));  // align Bytes
		}
		outs.write(reinterpret_cast<char*>(vals), sizeof(double)*nbMapSizes[i]);
	}
	outs.write(reinterpret_cast<char*>(nnbMapSizes),sizeof(int)*36000);
	for(int i=0; i<36000; i++) {
		int keys[nnbMapSizes[i]];
		double vals[nnbMapSizes[i]];
		int j = 0;
		for(auto& it1:nnbKeysEnergy[i]) {
			keys[j] = it1.first;
			vals[j] = it1.second;
			j++;
		}
		outs.write(reinterpret_cast<char*>(keys), sizeof(int)*nnbMapSizes[i]);
		if(nnbMapSizes[i]%2) {
			outs.write(reinterpret_cast<char*>(z0), sizeof(int));  // align Bytes
		}
		outs.write(reinterpret_cast<char*>(vals), sizeof(double)*nnbMapSizes[i]);
	}

	return EXIT_SUCCESS;
}

int BasePair6DEnergyTableCG::load() {
	ifstream ins;
	string fileName = datapath() + "../binaryCache/BasePair6DEnergyTableCG";
	ins.open(fileName, ios::in|ios::binary);
	if(!ins.is_open()) {
		throw("[Error]Fail to open " + fileName);
	}
	int tint[4];
	cm2Key.load(ins);
	ins.read(reinterpret_cast<char*>(&wtNb), sizeof(double));
	ins.read(reinterpret_cast<char*>(&wtNnb), sizeof(double));
	int nbMapSizes[36000], nnbMapSizes[36000];
	ins.read(reinterpret_cast<char*>(nbMapSizes),sizeof(int)*36000);
	for(int i=0; i<36000; i++) {
		int keys[nbMapSizes[i]];
		double vals[nbMapSizes[i]];
		if(nbMapSizes[i]%2) {
			int memSize = nbMapSizes[i]+1;
			int tint2[memSize];
			ins.read(reinterpret_cast<char*>(tint2), sizeof(int)*memSize);
			memcpy(keys, tint2, nbMapSizes[i]*sizeof(int));
		} else {
			ins.read(reinterpret_cast<char*>(keys), sizeof(int)*nbMapSizes[i]);
		}
		ins.read(reinterpret_cast<char*>(vals), sizeof(double)*nbMapSizes[i]);
		for(int j=0; j<nbMapSizes[i]; j++) {
			nbKeysEnergy[i].emplace(keys[j], vals[j]);
		}
	}
	ins.read(reinterpret_cast<char*>(nnbMapSizes),sizeof(int)*36000);
	for(int i=0; i<36000; i++) {
		int keys[nnbMapSizes[i]];
		double vals[nnbMapSizes[i]];
		if(nnbMapSizes[i]%2) {
			int memSize = nnbMapSizes[i]+1;
			int tint2[memSize];
			ins.read(reinterpret_cast<char*>(tint2), sizeof(int)*memSize);
			memcpy(keys, tint2, nnbMapSizes[i]*sizeof(int));
		} else {
			ins.read(reinterpret_cast<char*>(keys), sizeof(int)*nnbMapSizes[i]);
		}
		ins.read(reinterpret_cast<char*>(vals), sizeof(double)*nnbMapSizes[i]);
		for(int j=0; j<nnbMapSizes[i]; j++) {
			nnbKeysEnergy[i].emplace(keys[j], vals[j]);
		}
	}

	return EXIT_SUCCESS;
}

BasePair6DEnergyTableCG::~BasePair6DEnergyTableCG() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPpredNA */
