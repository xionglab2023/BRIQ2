/*
 * TestSingleSidechain.cpp
 *
 */

#include "protein/SingleSidechainModelingTemplate.h"
#include "para/ProParameter.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>

using namespace std;
using namespace NSPprotein;
using namespace NSPpara;
using namespace NSPtools;



double getRMS(const string datFile, ProEnergyTable* pe, ProParameter* para1, ProParameter* para2, ProParameter* para3, const string& scTag, double wtLocal, ofstream& out){

	ifstream file;
	file.open(datFile, ios::in);
	if(!file.is_open()){
		cout << "fail to open file: " << datFile << endl;
		exit(1);
	}

	ResName* rn = new ResName();
	cout << "init scLib" << endl;
	ResScRotamerLib* scLib = new ResScRotamerLib();
	cout << "init bbLib" << endl;
	ResBBRotamerLib* bbLib = new ResBBRotamerLib();

	AtomLib* atLib = new AtomLib();
	cout << "init ec" << endl;
	DesignPara* dp = new DesignPara("");
	EnergyCalculator* ec = new EnergyCalculator(pe, dp, para1);
	string s;
	vector<string> spt;
	vector<string> lines;
	vector<double> rmsList;
	int index = 1;
	char xx[200];
	while(getline(file, s)){
		if(s[0] == '#'){
			splitString(s, " ", &spt);
			int neighborNum = atoi(spt[2].c_str());
			lines.clear();
			for(int i=0;i<neighborNum+1;i++){
				getline(file, s);
				lines.push_back(s);
			}

			SingleSidechainModelingTemplate* st = new SingleSidechainModelingTemplate(lines, bbLib, scLib, rn);
			double rms = st->bestRotamerRMSD2(ec);
			sprintf(xx, "%d", index);

			string fileid = string(xx);
			index++;

			/*
			string out1 = "/user/xiongpeng/proteinModeling/singleSidechainModeling/out/p"+fileid+"-init.pdb";
			string out2 = "/user/xiongpeng/proteinModeling/singleSidechainModeling/out/p"+fileid+"-out.pdb";

			st->printPDB(out1, out2, atLib);
			*/

			out << rn->intToTri(st->targetNode->aaType) << " " << rms << endl;

			rmsList.push_back(rms);
			delete st;
		}
		if(rmsList.size() %10000 == 0){
			cout << rmsList.size() << endl;
		}
	}


	double meanRMS = 0;
	for(int i=0;i<rmsList.size();i++){
		meanRMS += rmsList[i];
	}
	meanRMS = meanRMS/rmsList.size();


	delete rn;
	delete scLib;
	delete bbLib;
	delete ec;

	sprintf(xx, "%8.6f", meanRMS);
	out << string(xx) << endl;
	return meanRMS;
}

void polarFinal(const string& polarPara1, const string& polarPara2, const string& polarPara3, const string& tri, const string& result){
	string datFile = "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/split/" + tri + ".dat";
	string filePolarS = "/user/xiongpeng/cpp/ProteinModeling/data/para/" + polarPara1;
	string filePolarM = "/user/xiongpeng/cpp/ProteinModeling/data/para/"+polarPara2;
	string filePolarL = "/user/xiongpeng/cpp/ProteinModeling/data/para/"+polarPara3;


	ProParameter* pa1 = new ProParameter("/user/xiongpeng/proteinModeling/para/qSR", filePolarS);
	ProParameter* pa2 = new ProParameter("/user/xiongpeng/proteinModeling/para/qMR", filePolarM);
	ProParameter* pa3 = new ProParameter("/user/xiongpeng/proteinModeling/para/qLR", filePolarL);
	string scTag = "20";
	ProEnergyTable* pe = new ProEnergyTable();

	ofstream out;
	out.open("/user/xiongpeng/proteinModeling/para/result/" + result, ios::out);
	double wtLocal = 1.0;
	getRMS(datFile, pe, pa1, pa2, pa3, scTag, wtLocal, out);
	out.close();
}

void vdwFinal(const string& vdw, const string& tri, const string& result){

	cout << "a" << endl;
	string datFile = "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/split/" + tri + ".dat";
	string filePolarS = "/user/xiongpeng/cpp/ProteinModeling/data/para/P0S";
	string filePolarM = "/user/xiongpeng/cpp/ProteinModeling/data/para/P0M";
	string filePolarL = "/user/xiongpeng/cpp/ProteinModeling/data/para/P0L";

	string fileVdw = "/user/xiongpeng/cpp/ProteinModeling/data/para/" + vdw;

	ProParameter* pa1 = new ProParameter(fileVdw, filePolarS);
	ProParameter* pa2 = new ProParameter(fileVdw, filePolarM);
	ProParameter* pa3 = new ProParameter(fileVdw, filePolarL);
	cout << "b" << endl;

	string scTag = "20";
	ProEnergyTable* pe = new ProEnergyTable();
	cout << "c" << endl;
	ofstream out;
	out.open("/user/xiongpeng/proteinModeling/para/result/" + result, ios::out);
	double wtLocal = 1.0;
	cout << "d" << endl;
	getRMS(datFile, pe, pa1, pa2, pa3, scTag, wtLocal, out);
	out.close();
}

void testPolarParaSep1(const string& initPara1, const string& initPara2, int polarPairID, double p1, double p2, double p3, const string& result) {
	map<int, string> datFileMap;
	ifstream file;
	string line;
	vector<string> spt;
	file.open("/user/xiongpeng/proteinModeling/singleSidechainModeling/data/polarData/pairToIndex", ios::in);
	while(getline(file, line)){
		splitString(line, " ", &spt);
		int pairID = atoi(spt[0].c_str());
		datFileMap[pairID] = "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/polarData/dat-" + spt[1];
	}

	if(datFileMap.find(polarPairID) == datFileMap.end()){
		cout << "invalid polarPairID: " << polarPairID << endl;
		exit(0);
	}

	string datFile = datFileMap[polarPairID];
	string filePolarS = "/user/xiongpeng/cpp/ProteinModeling/data/para/" + initPara1;
	string filePolarML = "/user/xiongpeng/cpp/ProteinModeling/data/para/"+initPara2;
	ProParameter* pa1 = new ProParameter("/user/xiongpeng/proteinModeling/para/qSR", filePolarS);
	ProParameter* pa2 = new ProParameter("/user/xiongpeng/proteinModeling/para/qMR", filePolarML);
	ProParameter* pa3 = new ProParameter("/user/xiongpeng/proteinModeling/para/qLR", filePolarML);


	pa1->paraPolar[0][polarPairID] = p1;
	pa1->paraPolar[1][polarPairID] = p2;
	pa1->paraPolar[2][polarPairID] = p3;

	/*
	pa2->paraPolar[0][polarPairID] = p1;
	pa2->paraPolar[1][polarPairID] = p2;
	pa2->paraPolar[2][polarPairID] = p3;

	pa3->paraPolar[0][polarPairID] = p1;
	pa3->paraPolar[1][polarPairID] = p2;
	pa3->paraPolar[2][polarPairID] = p3;
	*/

	string scTag = "20";

	ProEnergyTable* pe = new ProEnergyTable();

	ofstream out;
	out.open("/user/xiongpeng/proteinModeling/para/result/" + result, ios::out);
	double wtLocal = 1.0;
	getRMS(datFile, pe, pa1, pa2, pa3, scTag, wtLocal, out);
	out.close();
}

void testPolarParaSep23(const string& initPara1, const string& initPara2, int polarPairID, double p1, double p2, double p3, const string& result) {
	map<int, string> datFileMap;
	ifstream file;
	string line;
	vector<string> spt;
	file.open("/user/xiongpeng/proteinModeling/singleSidechainModeling/data/polarData/pairToIndex", ios::in);
	while(getline(file, line)){
		splitString(line, " ", &spt);
		int pairID = atoi(spt[0].c_str());
		datFileMap[pairID] = "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/polarData/dat-" + spt[1];
	}

	if(datFileMap.find(polarPairID) == datFileMap.end()){
		cout << "invalid polarPairID: " << polarPairID << endl;
		exit(0);
	}

	string datFile = datFileMap[polarPairID];
	string filePolarS = "/user/xiongpeng/cpp/ProteinModeling/data/para/" + initPara1;
	string filePolarML = "/user/xiongpeng/cpp/ProteinModeling/data/para/"+initPara2;
	ProParameter* pa1 = new ProParameter("/user/xiongpeng/proteinModeling/para/qSR", filePolarS);
	ProParameter* pa2 = new ProParameter("/user/xiongpeng/proteinModeling/para/qMR", filePolarML);
	ProParameter* pa3 = new ProParameter("/user/xiongpeng/proteinModeling/para/qLR", filePolarML);

	/*
	pa1->paraPolar[0][polarPairID] = p1;
	pa1->paraPolar[1][polarPairID] = p2;
	pa1->paraPolar[2][polarPairID] = p3;
	*/

	pa2->paraPolar[0][polarPairID] = p1;
	pa2->paraPolar[1][polarPairID] = p2;
	pa2->paraPolar[2][polarPairID] = p3;

	pa3->paraPolar[0][polarPairID] = p1;
	pa3->paraPolar[1][polarPairID] = p2;
	pa3->paraPolar[2][polarPairID] = p3;

	string scTag = "20";

	ProEnergyTable* pe = new ProEnergyTable();

	ofstream out;
	out.open("/user/xiongpeng/proteinModeling/para/result/" + result, ios::out);
	double wtLocal = 1.0;
	getRMS(datFile, pe, pa1, pa2, pa3, scTag, wtLocal, out);
	out.close();
}


void testDamping(const string& initDamp, int dampID, double value, const string& result){
	vector<string> pgList;
	pgList.push_back("BBO");
	pgList.push_back("BBN");
	pgList.push_back("GLN");
	pgList.push_back("GLU");
	pgList.push_back("ASP");
	pgList.push_back("ASN");
	pgList.push_back("SER");
	pgList.push_back("THR");
	pgList.push_back("LYS");
	pgList.push_back("ARG");
	pgList.push_back("HIS");
	pgList.push_back("PHE");
	pgList.push_back("TYR");
	pgList.push_back("TRP");

	int pgID = dampID%14;
	string datFile = "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/split/" + pgList[pgID] + ".dat";
	string filePolar = "/user/xiongpeng/cpp/ProteinModeling/data/para/p5";
	ProParameter* pa1 = new ProParameter("/user/xiongpeng/proteinModeling/para/qSR", filePolar);
	ProParameter* pa2 = new ProParameter("/user/xiongpeng/proteinModeling/para/qMR", filePolar);
	ProParameter* pa3 = new ProParameter("/user/xiongpeng/proteinModeling/para/qLR", filePolar);
	pa1->polarDamping[dampID] = value;
	pa2->polarDamping[dampID] = value;
	pa3->polarDamping[dampID] = value;

	string scTag = "20";

	ProEnergyTable* pe = new ProEnergyTable();

	ofstream out;
	out.open("/user/xiongpeng/proteinModeling/para/result/" + result, ios::out);
	double wtLocal = 1.0;
	getRMS(datFile, pe, pa1, pa2, pa3, scTag, wtLocal, out);
	out.close();
}

void test(const string& tri, const string& paraFile1, const string& paraFile2, const string& paraFile3, const string& scTag, double wtLocal, const string& result){

	cout << "init pa:" << endl;
	cout << paraFile1 << endl;
	ProParameter* pa1 = new ProParameter("/user/xiongpeng/proteinModeling/para/" + paraFile1);
	ProParameter* pa2 = new ProParameter("/user/xiongpeng/proteinModeling/para/" + paraFile2);
	ProParameter* pa3 = new ProParameter("/user/xiongpeng/proteinModeling/para/" + paraFile3);

	cout << "init pe:" << endl;
	ProEnergyTable* pe = new ProEnergyTable();

	ofstream out;
	out.open("/user/xiongpeng/proteinModeling/para/result/" + result, ios::out);

	string datFile = "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/split/" + tri + ".dat";
	getRMS(datFile, pe, pa1, pa2, pa3, scTag, wtLocal, out);
	out.close();
}


int main(int argc, char** argv){
	/*
	string tri = string(argv[1]);
	string paraFile1 = string(argv[2]);
	string paraFile2 = string(argv[3]);
	string paraFile3 = string(argv[4]);
	string scTag = string(argv[5]);
	double wtLocal = atof(argv[6]);
	string outFile = string(argv[7]);
	test(tri, paraFile1, paraFile2, paraFile3, scTag, wtLocal, outFile);
	*/

	/*
	string initPolar1 = string(argv[1]);
	string initPolar2 = string(argv[2]);
	int polarPairID = atoi(argv[3]);
	double p1 = atof(argv[4]);
	double p2 = atof(argv[5]);
	double p3 = atof(argv[6]);
	string result = string(argv[7]);
	int paraSep = atoi(argv[8]);

	if(paraSep == 1)
		testPolarParaSep1(initPolar1, initPolar2, polarPairID,  p1,  p2,  p3,  result);
	else if(paraSep == 2)
		testPolarParaSep23(initPolar1, initPolar2, polarPairID,  p1,  p2,  p3,  result);
	*/


	string vdw = string(argv[1]);
	string tri = string(argv[2]);
	string result = string(argv[3]);
	vdwFinal(vdw, tri, result);


	/*
	string initDamp = string(argv[1]);
	int dampID = atoi(argv[2]);
	double p = atof(argv[3]);
	string result = string(argv[4]);
	testDamping(initDamp, dampID, p, result);
	*/
}




