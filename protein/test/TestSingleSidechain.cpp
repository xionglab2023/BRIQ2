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



double getRMS(const string datFile, ProEnergyTable* pe, ProParameter* pa, const string& scTag, ofstream& out){

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
	EnergyCalculator* ec = new EnergyCalculator(pe, dp, pa);

	cout << "run:" << endl;
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

			//out << rn->intToTri(st->targetNode->aaType) << " " << rms << endl;

			//if(index > 50) break;

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

void polarFinal(const string& polarPara, const string& tri, const string& result){
	string datFile = "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/split/" + tri + ".dat";
	string filePolar = "/user/xiongpeng/cpp/ProteinModeling/data/para/"+polarPara;

	ProParameter* pa = new ProParameter("/user/xiongpeng/cpp/ProteinModeling/data/para/V0", filePolar);
	string scTag = "20";
	DesignPara* dp = new DesignPara("");

	ProEnergyTable* pe = new ProEnergyTable();

	ofstream out;
	out.open("/user/xiongpeng/proteinModeling/para/result/" + result, ios::out);
	getRMS(datFile, pe, pa, scTag, out);
	out.close();
}

void testVdwPara(const string& vdwPara, const string& polarPara, int vdwID, double value, const string& result) {

	ifstream file;
	string line;
	vector<string> spt;
	string datFile;

	char xx[200];

	/*
	if(vdwID < 5) {
		sprintf(xx, "%d", vdwID);
		datFile = "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/vdwData/wd/wd.dat";
	}
	else if(vdwID < 530) {
		int t = vdwID%5;
		sprintf(xx, "%d", (vdwID-5)/5);
		if(t < 2)
			datFile = "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/vdwData/S/dat-"+ string(xx);
		else
			datFile = "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/vdwData/L/dat-"+ string(xx);
	}
	else {
		cout << "invalid vdwID: "<< vdwID << endl;
		exit(0);
	}
*/



	if(vdwID < 5) {
		sprintf(xx, "%d", vdwID);
		datFile = "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/vdwData/wd/sub2/wd-"+string(xx)+".sub";
	}
	else if(vdwID < 530) {
		int t = vdwID%5;
		sprintf(xx, "%d", (vdwID-5)/5);
		if(t < 2)
			datFile = "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/vdwData/S/sub2/dat-"+ string(xx)+".sub";
		else
			datFile = "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/vdwData/L/sub2/dat-"+ string(xx)+".sub";
	}
	else {
		cout << "invalid vdwID: "<< vdwID << endl;
		exit(0);
	}



	string fileVdw = "/user/xiongpeng/cpp/ProteinModeling/data/para/"+vdwPara;
	string filePolar = "/user/xiongpeng/cpp/ProteinModeling/data/para/"+polarPara;

	cout << "init pa" << endl;
	ProParameter* pa = new ProParameter(fileVdw, filePolar);

	if(vdwID < 5) {
		pa->vdwSepWeight[vdwID] = value;
	}
	else {
		pa->paraVdw[vdwID-5] = value;
	}

	pa->updateVdwWellDepth();

	string scTag = "20";

	ProEnergyTable* pe = new ProEnergyTable();

	ofstream out;
	out.open("/user/xiongpeng/proteinModeling/para/result/" + result, ios::out);
	getRMS(datFile, pe, pa, scTag, out);
	out.close();
}


void testPolarPara(const string vdwPara,const string& polarPara, int paraIndex, double value, const string& result) {
	int polarPairID = paraIndex/6;
	int pid = paraIndex%6;
	string datFile;
	char xx[200];

	/*
	if(pid == 0 || pid == 1 || pid == 2)
	{
		sprintf(xx, "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/polarData/S/dat-%d", polarPairID);
		datFile = string(xx);
	}
	else if(pid == 3 || pid == 4 || pid == 5)
	{
		sprintf(xx, "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/polarData/M/dat-%d", polarPairID);
		datFile = string(xx);
	}
	else if(pid == 6 || pid == 7 || pid == 8)
	{
		sprintf(xx, "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/polarData/L/dat-%d", polarPairID);
		datFile = string(xx);
	}

	if(paraIndex >= 4914)
	{
		datFile = "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/polarData/damp/mix.dat";
	}
	*/


	if(pid == 0 || pid == 1 || pid == 2)
	{
		sprintf(xx, "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/polarData/S/sub2/dat-%d.sub", polarPairID);
		datFile = string(xx);
	}
	else if(pid == 3 || pid == 4 || pid == 5)
	{
		sprintf(xx, "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/polarData/L/sub2/dat-%d.sub", polarPairID);
		datFile = string(xx);
	}

	if(paraIndex >= 1428)
	{
		cout << "invalid paraID" << endl;
		exit(0);

		//sprintf(xx, "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/polarData/damp/sub/mix-%d.sub", paraIndex-4914);
		//datFile = string(xx);
	}


	string filePolar = "/user/xiongpeng/cpp/ProteinModeling/data/para/"+polarPara;
	string fileVdw = "/user/xiongpeng/cpp/ProteinModeling/data/para/" + vdwPara;
	ProParameter* pa = new ProParameter(fileVdw, filePolar);

	if(paraIndex < 1428)
		pa->paraPolar[polarPairID][pid] = value;
//	else if(paraIndex < 4942)
//		pa->polarDamping[paraIndex-4914] = value;
	else
	{
		cout << "invalid paraIndex: " << paraIndex << endl;
		exit(0);
	}

	string scTag = "20";

	ProEnergyTable* pe = new ProEnergyTable();

	ofstream out;
	out.open("/user/xiongpeng/proteinModeling/para/result/" + result, ios::out);
	getRMS(datFile, pe, pa, scTag, out);
	out.close();
}

void testDampPara(const string vdwPara,const string& polarPara, int dampID, double value, const string& result) {

	ifstream file;
	string line;
	vector<string> spt;
	string datFile;

	char xx[200];

	sprintf(xx, "%d", dampID);
	datFile = "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/polarData/damp/sub/mix-"+string(xx)+".sub";


	string fileVdw = "/user/xiongpeng/cpp/ProteinModeling/data/para/"+vdwPara;
	string filePolar = "/user/xiongpeng/cpp/ProteinModeling/data/para/"+polarPara;

	cout << "init pa" << endl;
	ProParameter* pa = new ProParameter(fileVdw, filePolar);

	pa->polarDamping[dampID] = value;
	pa->updateVdwWellDepth();

	string scTag = "20";

	ProEnergyTable* pe = new ProEnergyTable();

	ofstream out;
	out.open("/user/xiongpeng/proteinModeling/para/result/" + result, ios::out);
	getRMS(datFile, pe, pa, scTag, out);
	out.close();
}


int main(int argc, char** argv){

	string initVdw = string(argv[1]);
	string initPolar = string(argv[2]);

	int paraIndex = atoi(argv[3]);
	double paraValue = atof(argv[4]);

	string result = string(argv[5]);
	string tag = string(argv[6]);


	if(tag == "V")
		testVdwPara(initVdw, initPolar, paraIndex, paraValue, result);
	else if(tag == "P")
		testPolarPara(initVdw, initPolar, paraIndex,  paraValue, result);
	else if(tag == "D")
		testDampPara(initVdw, initPolar, paraIndex, paraValue, result);
}




