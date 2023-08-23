/*
 * TestScModeling.cpp
 *
 */


#include "protein/SidechainModelingTemplate.h"
#include "para/ProParameter.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>

using namespace std;
using namespace NSPprotein;
using namespace NSPpara;
using namespace NSPtools;


void predSidechain(){
	string fileVdw = "/user/xiongpeng/cpp/ProteinModeling/data/para/V4";
	ProParameter* pa1 = new ProParameter(fileVdw, "/user/xiongpeng/cpp/ProteinModeling/data/para/P4S");
	ProParameter* pa2 = new ProParameter(fileVdw, "/user/xiongpeng/cpp/ProteinModeling/data/para/P4M");
	ProParameter* pa3 = new ProParameter(fileVdw, "/user/xiongpeng/cpp/ProteinModeling/data/para/P4L");

	pa1->updateVdwWellDepth();
	pa2->updateVdwWellDepth();
	pa3->updateVdwWellDepth();

	string scTag = "20";

	ProEnergyTable* pe = new ProEnergyTable();
	DesignPara* dp = new DesignPara("");
	EnergyCalculator* ec = new EnergyCalculator(pe, dp, pa1);
	ifstream file;
	file.open("/user/xiongpeng/proteinModeling/sidechainModeling/DLPacker");
	string line;
	vector<string> pdbList;
	while(getline(file, line)){
		pdbList.push_back(line);
	}
	file.close();
	for(int i=0;i<pdbList.size();i++){
		string pdbFile = "/user/xiongpeng/proteinModeling/sidechainModeling/DLPacker/rawPDB/" + pdbList[i] + ".pdb";
		SidechainModelingTemplate* st = new SidechainModelingTemplate(pdbFile);
		st->runMC(ec);
		delete st;
	}
	delete ec;
	delete pe;
	delete pa1;
	delete pa2;
	delete pa3;
}

void testList(const string& pdbListFile, const string& outFile, double dShiftS, double dLamdaS, double dShiftL, double dLamdaL){
	ifstream f;
	f.open(pdbListFile.c_str(), ios::in);
	string line;
	string fileVdw = "/user/xiongpeng/cpp/ProteinModeling/data/para/V4";
	string filePolar = "/user/xiongpeng/cpp/ProteinModeling/data/para/P4";

	ProParameter* pa = new ProParameter(fileVdw, filePolar);
	for(int i=0;i<105;i++){
		pa->paraVdw[i*5] += dShiftS;
		pa->paraVdw[i*5+1] += dLamdaS;
		pa->paraVdw[i*5+2] += dShiftL;
		pa->paraVdw[i*5+3] += dLamdaL;
	}
	pa->updateVdwWellDepth();

	string scTag = "20";

	ProEnergyTable* pe = new ProEnergyTable();
	DesignPara* dp = new DesignPara("");
	EnergyCalculator* ec = new EnergyCalculator(pe, dp, pa);
	ofstream out;
	out.open(outFile.c_str(), ios::out);

	while(getline(f, line)){
		cout << line << endl;
		string pdbFile = "/user/xiongpeng/lib/protein/pdbs/"+line+".pdb";
		SidechainModelingTemplate* st = new SidechainModelingTemplate(pdbFile);
		st->runMC(ec);
		st->printRMS(out);
		delete st;
	}
	out.close();
	f.close();
}

void testPDB(const string& pdbFile, const string& outFile, double dShiftS, double dLamdaS, double dShiftL, double dLamdaL){

	cout << "start" << endl;
	SidechainModelingTemplate* st = new SidechainModelingTemplate(pdbFile);

	string fileVdw = "/user/xiongpeng/cpp/ProteinModeling/data/para/V4";
	string filePolar = "/user/xiongpeng/cpp/ProteinModeling/data/para/P4";

	ProParameter* pa = new ProParameter(fileVdw, filePolar);


	for(int i=0;i<105;i++){
		pa->paraVdw[i*5] += dShiftS;
		pa->paraVdw[i*5+1] += dLamdaS;
		pa->paraVdw[i*5+2] += dShiftL;
		pa->paraVdw[i*5+3] += dLamdaL;
		printf("%-6.3f %5.3f %6.3f %5.3f %6.3f\n", pa->paraVdw[i*5],pa->paraVdw[i*5+1],pa->paraVdw[i*5+2],pa->paraVdw[i*5+3],pa->paraVdw[i*5+4]);
	}


	pa->updateVdwWellDepth();
	string scTag = "20";

	ProEnergyTable* pe = new ProEnergyTable();
	DesignPara* dp = new DesignPara("");
	EnergyCalculator* ec = new EnergyCalculator(pe, dp, pa);

	cout << "runMC" << endl;
	st->runMC(ec);
	cout << "finish" << endl;
	st->printPDB(outFile);
	delete st;
}


int main(int argc, char** argv){


	string pdbList = string(argv[1]);
	string outFile = string(argv[2]);
	double dShiftS = atof(argv[3]);
	double dLamdaS = atof(argv[4]);
	double dShiftL = atof(argv[5]);
	double dLamdaL = atof(argv[6]);
	testList(pdbList, outFile, dShiftS, dLamdaS, dShiftL, dLamdaL);

	/*
	string pdbFile = string(argv[1]);
	string outFile = string(argv[2]);
	double dShiftS = atof(argv[3]);
	double dLamdaS = atof(argv[4]);
	double dShiftL = atof(argv[5]);
	double dLamdaL = atof(argv[6]);
	testPDB(pdbFile, outFile, dShiftS, dLamdaS, dShiftL, dLamdaL);
	*/
}

