/*
 * RiboseOxygenEnergyTable.cpp
 *
 */

#include "forcefield/RiboseOxygenEnergyTable.h"

namespace NSPforcefield {

RiboseOxygenEnergyTable::RiboseOxygenEnergyTable() {
	// TODO Auto-generated constructor stub
	string path = NSPdataio::datapath() + "rna/riboseOxygen/";
	vector<string> oxygenNames;

	oxygenNames.push_back("O3");
	oxygenNames.push_back("O4");
	oxygenNames.push_back("O5");
	string augc = "AUGCatgc";

	vector<string> spt;
	double ene;
	int directIndex;
	for(int i=0;i<8;i++){
		for(int j=0;j<3;j++){
			string energyFile = path + augc.substr(i,1) + "-" + oxygenNames[j]+".ene";
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
				this->etList[i*3+j].push_back(ene);
			}
			file.close();
		}
	}

	/*
	 * sep = -1
	 */

	vector<string> oxygenNamesM1;
	oxygenNamesM1.push_back("O4");
	for(int i=0;i<8;i++){
		string energyFile = path + "sep-1/" + augc.substr(i,1) + "-" + oxygenNamesM1[0]+".ene";
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
			this->etListM1[i].push_back(ene);
		}
		file.close();
	}

}



RiboseOxygenEnergyTable::~RiboseOxygenEnergyTable() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPforcefield */
