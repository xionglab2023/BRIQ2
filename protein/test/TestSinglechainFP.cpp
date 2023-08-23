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
			//st->printNativeEnergy(ec);
			double rms = st->bestRotamerRMSD2(ec);
			sprintf(xx, "%d", index);

			string fileid = string(xx);
			index++;


			/*
			string out1 = "/user/xiongpeng/proteinModeling/singleSidechainModeling/out/p"+fileid+"-init.pdb";
			string out2 = "/user/xiongpeng/proteinModeling/singleSidechainModeling/out/p"+fileid+"-out.pdb";
			st->printPDB(out1, out2, atLib);
			cout << "model: " << index << endl;
			st->printNativeEnergy(ec);
			cout << rn->intToTri(st->targetNode->aaType) << " " << rms << endl;
			if(index > 20) break;
			*/


			out << rn->intToTri(st->targetNode->aaType) << " " << rms << endl;
			//if(index == 2) break;



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

void polarFinal(const string& polarPara , const string& vdwPara, const string& tri, const string& result){
	string datFile = "/user/xiongpeng/proteinModeling/singleSidechainModeling/data/split/" + tri + ".dat";
	string filePolar = "/user/xiongpeng/cpp/ProteinModeling/data/para/" + polarPara;
	string fileVdw = "/user/xiongpeng/cpp/ProteinModeling/data/para/"+ vdwPara;

	ProParameter* pa = new ProParameter(fileVdw, filePolar);
	string scTag = "20";
	ProEnergyTable* pe = new ProEnergyTable();

	ofstream out;
	out.open("/user/xiongpeng/proteinModeling/para/result/" + result, ios::out);
	getRMS(datFile, pe, pa, scTag, out);
	out.close();
}


int main(int argc, char** argv){

	string polar = string(argv[1]);
	string vdw = string(argv[2]);
	string tri = string(argv[3]);
	string result = string(argv[4]);

	polarFinal(polar, vdw,  tri, result);

}




