/*
 * TestDesign.cpp
 *
 *  Created on: 2022Äê1ÔÂ9ÈÕ
 *      Author: pengx
 */

#include "protein/SeqDesignTemplate.h"
#include "protein/SeqDesignTemplateFast.h"
#include "para/ProParameter.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>

using namespace std;
using namespace NSPprotein;
using namespace NSPpara;
using namespace NSPtools;


void scRotamerClustering() {
	ResScRotamerLib* scLib = new ResScRotamerLib();


	ResName rn;
	for(int i=0;i<20;i++){
		vector<ResScRotamer*> rotList = scLib->rotList[i];

		vector<ResScRotamer*> clusterCenters;
		vector<vector<int>> members;

		int clusterNum = scLib->getAminoAcidRotamerClusterNum(i);

		if(clusterNum == 0) continue;
		ofstream out;
		out.open("/user/xiongpeng/cpp/ProteinModeling/data/scRotamer/rot/twoLevelUnique/"+rn.intToTri(i)+".cluster", ios::out);

		for(int j=0;j<clusterNum;j++){
			int rotIndex = scLib->rotCluster[i][j][1];
			ResScRotamer* rot = rotList[rotIndex];
			clusterCenters.push_back(rot);
			members.push_back(vector<int>());
			members[j].push_back(rotIndex);
		}

		double d;
		for(int j=0;j<rotList.size();j++){
			double minD = 999.9;
			int minIndex = -1;
			for(int k=0;k<clusterNum;k++){
				d = rotList[j]->distanceTo(clusterCenters[k]);
				if(d < minD){
					minD = d;
					minIndex = k;
				}
			}
			if(j != members[minIndex][0])
				members[minIndex].push_back(j);
		}

		for(int j=0;j<clusterNum;j++){
			vector<int> mem = members[j];
			for(int k=0;k<mem.size();k++){
				out << mem[k] << " ";
			}
			out << endl;
		}
		out.close();

	}
}

void printS2(const string& pdbID) {
	clock_t start = clock();
	srand(time(0));
	ProParameter* pa = new ProParameter();
	DesignPara* dp = new DesignPara("/user/xiongpeng/sword/designPara/py");
	ProEnergyTable* pe = new ProEnergyTable();
	ProS1S2Energy* eS1S2= new ProS1S2Energy(dp, 0);
	EnergyCalculator* ec = new EnergyCalculator(pe, dp, pa);

	string inputFile = "/user/xiongpeng/sword/train/" + pdbID + "/input";

	string s2out = "/user/xiongpeng/sword/train/S2/"+pdbID+".s2";

	ResBBRotamerLib* bbLib = new ResBBRotamerLib();
	ResScRotamerLibMini* scLib = new ResScRotamerLibMini();
	AtomLib* atLib = new AtomLib();
	SeqDesignTemplateFast* dt = new SeqDesignTemplateFast(inputFile, bbLib, scLib, atLib, ec, eS1S2);
	dt->printS2(s2out);
	delete dt;
}

void train(const string& pdbID, const string& id, const string& bpIndex){

	//string filePolar = "/user/xiongpeng/cpp/ProteinModeling/data/para/pDes";
	//string designPara = "/user/xiongpeng/cpp/ProteinModeling/data/para/design/" + paraFile;

	clock_t start = clock();
	srand(time(0));

	cout << "train: " << pdbID << endl;
	ProParameter* pa = new ProParameter();
	DesignPara* dp = new DesignPara("/user/xiongpeng/sword/designPara/bp/p" + id);
	cout << "PE" << endl;
	ProEnergyTable* pe = new ProEnergyTable();
	cout << "s1s2" << endl;
	ProS1S2Energy* eS1S2= new ProS1S2Energy(dp, 0);
	cout << "ec" << endl;
	EnergyCalculator* ec = new EnergyCalculator(pe, dp, pa);

	string inputFile = "/user/xiongpeng/sword/train/" + pdbID + "/input";
	string output = "/user/xiongpeng/sword/train/designSeq"+ bpIndex + "/"+pdbID+".seq";

	ofstream of;
	of.open(output.c_str(), ios::out);
	ResBBRotamerLib* bbLib = new ResBBRotamerLib();
	ResScRotamerLibMini* scLib = new ResScRotamerLibMini();
	AtomLib* atLib = new AtomLib();

	string s2File = "/user/xiongpeng/sword/train/S2/"+pdbID+".s2";
	SeqDesignTemplateFast* dt = new SeqDesignTemplateFast(inputFile, bbLib, scLib, atLib, ec, eS1S2, s2File);
	dt->loadEAEM();

	clock_t t1 = clock();
	of << "time: " << (float)(t1-start)/CLOCKS_PER_SEC << "s" << endl;

	char xx[200];
	for(int i=0;i<5;i++) {

		dt->designMC();
		of << dt->getDesignSequence() << endl;
		sprintf(xx, "/user/xiongpeng/sword/train/designPDB%s/%s-%d.pdb",bpIndex.c_str(), pdbID.c_str(), i);
		string outpdb = string(xx);
		dt->printDesignPDB(outpdb);
	}
	of.close();

}

void testDes() {
	string filePolar = "/user/xiongpeng/cpp/ProteinModeling/data/para/pDes";
	string fileVdw = "/user/xiongpeng/cpp/ProteinModeling/data/para/vDes";

	srand(time(0));
	ProParameter* pa = new ProParameter(fileVdw, filePolar);
	pa->polarSepWeight[0] = 0.2;
	pa->polarSepWeight[1] = 0.8;
	pa->polarSepWeight[2] = 1.1;
	pa->polarSepWeight[3] = 1.1;
	pa->polarSepWeight[4] = 1.2;

	pa->wtRot = 1.0;


	pa->updateVdwWellDepth();

	clock_t start = clock();

	DesignPara* dp = new DesignPara("");
	ProEnergyTable* pe = new ProEnergyTable(dp);
	ProS1S2Energy* eS1S2= new ProS1S2Energy(dp, 1);
	EnergyCalculator* ec = new EnergyCalculator(pe, dp, pa);

	string inputFile = "/user/xiongpeng/proteinModeling/design/input";
	string output = "/user/xiongpeng/proteinModeling/design/des.pdb";

	ResBBRotamerLib* bbLib = new ResBBRotamerLib();
	ResScRotamerLib* scLib = new ResScRotamerLib();
	AtomLib* atLib = new AtomLib();
	SeqDesignTemplate* dt = new SeqDesignTemplate(inputFile,bbLib, scLib, atLib, ec, eS1S2);
	//dt->printResInfo();
	//dt->getPositionRanks();
	cout << "run" << endl;
	dt->runMC();
	//dt->printPDB(output);
	//dt->printDetailEnergy();
	//dt->printInvolvedPairInfo();
}

void design1ubq(){
	clock_t start = clock();
	srand(time(0));
	ProParameter* pa = new ProParameter();
	string paraFile = "/user/xiongpeng/sword/designPara/tp13";
	DesignPara* dp = new DesignPara(paraFile);
	ProEnergyTable* pe = new ProEnergyTable();
	ProS1S2Energy* eS1S2= new ProS1S2Energy(dp, 0);

	EnergyCalculator* ec = new EnergyCalculator(pe, dp, pa);

	string inputFile = "/user/xiongpeng/input";
	string output = "/user/xiongpeng/1ubq-desSeq";
	string outpdb = "/user/xiongpeng/1ubq.des.pdb";
	ofstream of;
	of.open(output.c_str(), ios::out);

	ResBBRotamerLib* bbLib = new ResBBRotamerLib();
	ResScRotamerLibMini* scLib = new ResScRotamerLibMini();
	AtomLib* atLib = new AtomLib();
	SeqDesignTemplateFast* dt = new SeqDesignTemplateFast(inputFile, bbLib, scLib, atLib, ec, eS1S2);
	dt->loadEAEM();

	clock_t t1 = clock();
	cout << "time: " << (float)(t1-start)/CLOCKS_PER_SEC << "s" << endl;

	dt->printNatEnergy();

	dt->designMC();
	clock_t t2 = clock();
	cout << "time: "<< "des " << (float)(t2-start)/CLOCKS_PER_SEC << "s" << endl;
	of << dt->getDesignSequence() << endl;
	dt->printDesignPDB(outpdb);
	of.close();

	delete bbLib;
	delete scLib;
	delete atLib;
	delete dt;
}

int main(int argc, char** argv){

	/*
	string pdbID = string(argv[1]);
	printS2(pdbID);
	*/


	string pdbID = string(argv[1]);
	string id = string(argv[2]); //0
	string bpIndex = string(argv[3]);
	train(pdbID, id, bpIndex);


	//design1ubq();

	//string paraFile = string(argv[2]);
	/*
	float wtSurf = atof(argv[2]); //0.05
	float wtCore = atof(argv[3]); //0.5
	float wtPolar = atof(argv[4]); //1.0
	float wtRot = atof(argv[5]); //1.0
	float wtS1 = atof(argv[6]);  //4.0
	float wtS2 = atof(argv[7]); //2.2
	*/

}
