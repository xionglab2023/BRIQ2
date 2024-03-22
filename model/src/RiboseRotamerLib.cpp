/*
 * RiboseRotamerLib.cpp
 *
 *  Created on: 2022��9��6��
 *      Author: pengx
 */

#include "model/RiboseRotamerLib.h"

namespace NSPmodel {

RiboseRotamerLib::RiboseRotamerLib() {
	// TODO Auto-generated constructor stub
	string augc = "AUGCatgc";
	string path = NSPdataio::datapath();

	ifstream file;
	string s;
	int lv1 = 0;

	for(int i=0;i<8;i++){
		string baseType = augc.substr(i,1);
		string fileName = path+"riboseRotamer/"+baseType+".rot";
		file.open(fileName.c_str(), ios::in);
		if(! file.is_open()) {
			cout << "fail to open rotamer library file " << fileName << endl;
			exit(1);
		}
		int index = 0;

		while(getline(file,s)){
			RiboseRotamer* rot = new RiboseRotamer(s);
			rot->resType = i;

			if(i<4) { 
				//RNA
				if(index < 600)
					lv1 = 0;
				else if(index < 900)
					lv1 = 1;
				else if(index < 1200)
					lv1 = 2;
				else 
					lv1 = 3;
			}
			else {
				//DNA
				if(index < 300)
					lv1 = 0;
				else if(index < 375)
					lv1 = 1;
				else if(index < 1425)
					lv1 = 2;
				else
					lv1 = 3;
			}

			rot->rotTypeLv1 = lv1;
			rot->rotType = index;
			rotLib[i][index] = rot;
			rotLibCG[i][index] = new RiboseRotamerCG(rot);
			index++;
		}
		file.close();
	}

	for(int i=0;i<8;i++){
		string baseType = augc.substr(i,1);
		string fileName = path+"riboseRotamer/rot100/"+baseType+".rot";
		file.open(fileName.c_str(), ios::in);
		if(! file.is_open()) {
			cout << "fail to open rotamer library file " << fileName << endl;
			exit(1);
		}
		int index = 0;
		while(getline(file,s)){
			RiboseRotamer* rot = new RiboseRotamer(s);
			if(i<4) { 
				//RNA
				if(index < 40)
					lv1 = 0;
				else if(index < 60)
					lv1 = 1;
				else if(index < 80)
					lv1 = 2;
				else 
					lv1 = 3;
			}
			else {
				//DNA
				if(index < 20)
					lv1 = 0;
				else if(index < 25)
					lv1 = 1;
				else if(index < 95)
					lv1 = 2;
				else
					lv1 = 3;
			}
			rot->resType = i;
			rot->rotTypeLv1 = lv1;
			rot->rotType = index;
			rotLib100[i][index] = rot;
			rotLibCG100[i][index] = new RiboseRotamerCG(rot);
			index++;
		}
		file.close();
	}

	for(int i=0;i<8;i++){
		string baseType = augc.substr(i,1);
		string fileName = path+"riboseRotamer/rot20/"+baseType+".rot";
		file.open(fileName.c_str(), ios::in);
		if(! file.is_open()) {
			cout << "fail to open rotamer library file " << fileName << endl;
			exit(1);
		}
		int index = 0;
		while(getline(file,s)){
			RiboseRotamer* rot = new RiboseRotamer(s);
			if(i<4) { 
				//RNA
				if(index < 8)
					lv1 = 0;
				else if(index < 12)
					lv1 = 1;
				else if(index < 16)
					lv1 = 2;
				else 
					lv1 = 3;
			}
			else {
				//DNA
				if(index < 4)
					lv1 = 0;
				else if(index < 5)
					lv1 = 1;
				else if(index < 19)
					lv1 = 2;
				else
					lv1 = 3;
			}
			rot->resType = i;
			rot->rotTypeLv1 = lv1;
			rot->rotType = index;
			rotLib20[i][index] = rot;
			rotLibCG20[i][index] = new RiboseRotamerCG(rot);
			index++;
		}
		file.close();
	}
}

RiboseRotamerLib::RiboseRotamerLib(ForceFieldPara* para){
	string augc = "AUGCatgc";
	string path = NSPdataio::datapath();
	ifstream file;
	string s;
	int lv1;
	for(int i=0;i<8;i++){
		string baseType = augc.substr(i,1);
		string fileName = path+"riboseRotamer/"+baseType+".rot";
		file.open(fileName.c_str(), ios::in);
		if(! file.is_open()) {
			cout << "fail to open rotamer library file " << fileName << endl;
			exit(1);
		}

		int aaType = i%4;
		int index = 0;
		while(getline(file,s)){
			RiboseRotamer* rot = new RiboseRotamer(s);

			if(i < 4){
				if(index < 600)
					rot->energy = (rot->energy+para->rnaRiboseRotamerShift[aaType*4 + 0])*para->wtRibose;
				else if(index < 900)
					rot->energy = (rot->energy+para->rnaRiboseRotamerShift[aaType*4 + 1])*para->wtRibose;
				else if(index < 1200)
					rot->energy = (rot->energy+para->rnaRiboseRotamerShift[aaType*4 + 2])*para->wtRibose;
				else
					rot->energy = (rot->energy+para->rnaRiboseRotamerShift[aaType*4 + 3])*para->wtRibose;
			}

			if(i<4) { 
				//RNA
				if(index < 600)
					lv1 = 0;
				else if(index < 900)
					lv1 = 1;
				else if(index < 1200)
					lv1 = 2;
				else 
					lv1 = 3;
			}
			else {
				//DNA
				if(index < 300)
					lv1 = 0;
				else if(index < 375)
					lv1 = 1;
				else if(index < 1425)
					lv1 = 2;
				else
					lv1 = 3;
			}

			rot->resType = i;
			rot->rotTypeLv1 = lv1;
			rot->rotType = index;
			rotLib[i][index] = rot;
			rotLibCG[i][index] = new RiboseRotamerCG(rot);
			index++;
		}
		file.close();

	}

	for(int i=0;i<8;i++){
		string baseType = augc.substr(i,1);
		string fileName = path+"riboseRotamer/rot100/"+baseType+".rot";
		file.open(fileName.c_str(), ios::in);
		if(! file.is_open()) {
			cout << "fail to open rotamer library file " << fileName << endl;
			exit(1);
		}
		int aaType = i%4;
		int index = 0;
		while(getline(file,s)){
			RiboseRotamer* rot = new RiboseRotamer(s);
			if(i < 4){
				if(index < 40)
					rot->energy = (rot->energy+para->rnaRiboseRotamerShift[aaType*4 + 0])*para->wtRibose;
				else if(index < 60)
					rot->energy = (rot->energy+para->rnaRiboseRotamerShift[aaType*4 + 1])*para->wtRibose;
				else if(index < 80)
					rot->energy = (rot->energy+para->rnaRiboseRotamerShift[aaType*4 + 2])*para->wtRibose;
				else
					rot->energy = (rot->energy+para->rnaRiboseRotamerShift[aaType*4 + 3])*para->wtRibose;
			}

			if(i<4) { 
				//RNA
				if(index < 40)
					lv1 = 0;
				else if(index < 60)
					lv1 = 1;
				else if(index < 80)
					lv1 = 2;
				else 
					lv1 = 3;
			}
			else {
				//DNA
				if(index < 20)
					lv1 = 0;
				else if(index < 25)
					lv1 = 1;
				else if(index < 95)
					lv1 = 2;
				else
					lv1 = 3;
			}
			rot->resType = i;
			rot->rotTypeLv1 = lv1;
			rot->rotType = index;
			rotLib100[i][index] = rot;
			rotLibCG100[i][index] = new RiboseRotamerCG(rot);
			index++;
		}
		file.close();
	}

	for(int i=0;i<8;i++){
		string baseType = augc.substr(i,1);
		string fileName = path+"riboseRotamer/rot20/"+baseType+".rot";
		file.open(fileName.c_str(), ios::in);
		if(! file.is_open()) {
			cout << "fail to open rotamer library file " << fileName << endl;
			exit(1);
		}
		int index = 0;
		while(getline(file,s)){
			RiboseRotamer* rot = new RiboseRotamer(s);
			if(i<4) { 
				//RNA
				if(index < 8)
					lv1 = 0;
				else if(index < 12)
					lv1 = 1;
				else if(index < 16)
					lv1 = 2;
				else 
					lv1 = 3;
			}
			else {
				//DNA
				if(index < 4)
					lv1 = 0;
				else if(index < 5)
					lv1 = 1;
				else if(index < 19)
					lv1 = 2;
				else
					lv1 = 3;
			}
			rot->resType = i;
			rot->rotTypeLv1 = lv1;
			rot->rotType = index;
			rotLib20[i][index] = rot;
			rotLibCG20[i][index] = new RiboseRotamerCG(rot);
			index++;
		}
		file.close();
	}
}

RiboseRotamerLib::~RiboseRotamerLib() {
	// TODO Auto-generated destructor stub
	for(int i=0;i<8;i++){
		for(int j=0;j<1500;j++){
			delete rotLib[i][j];
			delete rotLibCG[i][j];
		}
	}

	for(int i=0;i<8;i++){
		for(int j=0;j<100;j++){
			delete rotLib100[i][j];
			delete rotLibCG100[i][j];
		}
	}

	for(int i=0;i<8;i++){
		for(int j=0;j<20;j++){
			delete rotLib20[i][j];
			delete rotLibCG20[i][j];
		}
	}	
}

} /* namespace NSPmodel */