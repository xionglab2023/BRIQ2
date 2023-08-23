/*
 * ProteinRNAEnergyTable.cpp
 *
 *  Created on: 2022Äê7ÔÂ19ÈÕ
 *      Author: pengx
 */

#include "forcefield/ProteinRNAEnergyTable.h"

namespace NSPforcefield {

ProteinRNAEnergyTable::ProteinRNAEnergyTable() {
	// TODO Auto-generated constructor stub
	string path = NSPdataio::datapath() + "complexPR/basePolarEne/";
	vector<string> polarNames;
	polarNames.push_back("BBN");
	polarNames.push_back("BBO");

	this->atomLib = new AtomLib();

	string augc = "AUGC";
	this->para = new XPara();
	vector<string> spt;
	double ene;

	int directIndex;
	for(int i=0;i<4;i++){
		for(int j=0;j<2;j++){
			string energyFile = path + augc.substr(i,1) + "-" + polarNames[j]+".ene";
			ifstream file;
			file.open(energyFile, ios::in);
			if(file.fail()) {
				cout << "can't open file: " << energyFile << endl;
				exit(0);
			}
			string line;

			while(getline(file, line)){
				NSPtools::splitString(line, " ", &spt);
				ene = atof(spt[0].c_str());
				ene = energyRescale(ene);
				this->etList[i*2+j].push_back(ene);
				if(spt.size() < 2) {
					this->directDep[i*2+j].push_back(false);
					XYZ t;
					this->hbAtom[i*2+j].push_back(t);
					this->ang1List[i*2+j].push_back(0.0);
					this->ang2List[i*2+j].push_back(0.0);
					this->ang3List[i*2+j].push_back(0.0);
				}
				else {
					this->directDep[i*2+j].push_back(true);
					XYZ t(atof(spt[1].c_str()), atof(spt[2].c_str()), atof(spt[3].c_str()));
					double ang1 = atof(spt[4].c_str());
					double ang2 = atof(spt[5].c_str());
					double ang3 = atof(spt[6].c_str());
					this->hbAtom[i*2+j].push_back(t);
					this->ang1List[i*2+j].push_back(ang1);
					this->ang2List[i*2+j].push_back(ang2);
					this->ang3List[i*2+j].push_back(ang3);
				}
			}
			file.close();
		}
	}
}

ProteinRNAEnergyTable::~ProteinRNAEnergyTable() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPforcefield */
