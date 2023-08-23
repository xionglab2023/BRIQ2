/*
 * ResPeptideRotamerLib.cpp
 *
 */

#include "model/ResPeptideRotamerLib.h"

namespace NSPmodel {
ResPeptideRotamerLib::ResPeptideRotamerLib() {
	string path = NSPdataio::datapath();
	ResName rn;

	ifstream file;
	string s;
	int index;
	vector<string> spt;
	float e;

	for(int i=0;i<20;i++){
		string triA = rn.intToTri(i);
		for(int j=0;j<20;j++){
			string triB = rn.intToTri(j);

			string libFile = path + "/pepRotamer/rots/" + triA + "-" + triB + ".rot";
			file.open(libFile.c_str(), ios::in);
			if(! file.is_open()) {
				cout << "fail to open peptide rotamer library file " << libFile << endl;
				exit(1);
			}
			int index = 0;
			getline(file, s);
			while(getline(file,s)){
				this->pepLib[i*1600+j*80 + index] = new ResPeptideRotamer(s);
				this->pepLib[i*1600+j*80 + index]->setRotID(index);
				index++;
			}
			file.close();


			libFile = path + "/pepRotamer/ene/" + triA + "-" + triB + ".ene";
			file.open(libFile.c_str(), ios::in);
			if(! file.is_open()) {
				cout << "fail to open peptide energy library file " << libFile << endl;
				exit(1);
			}

			index = 0;
			while(getline(file, s)){
				splitString(s, " ", &spt);

				for(int k=0;k<80;k++){
					e = atoi(spt[k].c_str());
					this->ene[(i*20+j)*80 + k][index] = e;
				}
				index ++;
			}
			file.close();


			libFile = path + "/pepRotamer/psiphi.lib";
			file.open(libFile.c_str(), ios::in);
			if(! file.is_open()) {
				cout << "fail to open psiphi library file " << libFile << endl;
				exit(1);
			}
			index = 0;
			while(getline(file, s)){
				splitString(s, " ", &spt);
				ppLibPhi[index] = atof(spt[0].c_str());
				ppLibPsi[index] = atof(spt[1].c_str());
				index++;
			}
			file.close();

			libFile = path + "/pepRotamer/psiphi.index";
			file.open(libFile.c_str(), ios::in);
			if(! file.is_open()) {
				cout << "fail to open psiphi index file " << libFile << endl;
				exit(1);
			}

			index = 0;
			while(getline(file, s)){
				psiphiToPPIndex[index] = atoi(s.c_str());
				index++;
			}
			file.close();

		}
	}

}

ResPeptideRotamerLib::~ResPeptideRotamerLib() {
	// TODO Auto-generated destructor stub
	for(int i=0;i<32000;i++){
		delete pepLib[i];
	}
}

}
