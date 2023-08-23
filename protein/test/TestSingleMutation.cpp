
#include "protein/BuildMutationTemplate.h"
#include "protein/SeqDesignTemplate.h"
#include "para/ProParameter.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>

using namespace std;
using namespace NSPprotein;
using namespace NSPpara;
using namespace NSPtools;

void readS1S2(const string& datFile, const string& outFile) {
	ifstream file;
	ofstream out;
	out.open(outFile, ios::out);
	DesignPara* dp = new DesignPara("");
	ProS1S2Energy* eS1S2= new ProS1S2Energy(dp, 0);

	file.open(datFile, ios::in);
	if(!file.is_open()){
		cout << "fail to open file: " << datFile << endl;
		exit(1);
	}

	ResName* rn = new ResName();
	cout << "init scLib" << endl;
	ResScRotamerLib* scLib = new ResScRotamerLib("");
	ResScRotamerLibMini* scLibSimp = new ResScRotamerLibMini();
	cout << "init bbLib" << endl;
	ResBBRotamerLib* bbLib = new ResBBRotamerLib();
	AtomLib* atLib = new AtomLib();

	string s;
	vector<string> spt;
	vector<string> lines;
	vector<int> rankList;
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

			BuildMutationTemplate* mt = new BuildMutationTemplate(lines, bbLib, scLib, scLibSimp, atLib, rn);
			mt->loadS1S2(eS1S2);
			mt->printS1S2(out);
			sprintf(xx, "%d", index);

			string fileid = string(xx);
			index++;
			delete mt;
		}
		if(index %1000 == 0){
			cout << index << endl;
		}
	}
	out.close();
	file.close();
}

void testSingleDesign(const string& datFile, ProEnergyTable* pe, ProS1S2Energy* eS1S2, DesignPara* dp, ofstream& out, float wtSurf, float wtCore, float wtPolar, double* result){
	//result: [0~19] amino acid count
	//result: [20] total count
	//result: [21] identital count
	//result: [22~25] deltaVol

	ifstream file;
	file.open(datFile, ios::in);
	if(!file.is_open()){
		cout << "fail to open file: " << datFile << endl;
		exit(1);
	}

	ProParameter* pa = new ProParameter();

	ResName* rn = new ResName();
	ResScRotamerLib* scLib = new ResScRotamerLib("");
	ResScRotamerLibMini* scLibSimp = new ResScRotamerLibMini();
	ResBBRotamerLib* bbLib = new ResBBRotamerLib();
	AtomLib* atLib = new AtomLib();
	EnergyCalculator* ec = new EnergyCalculator(pe, dp, pa);

	string s;
	vector<string> spt;
	vector<string> lines;
	vector<int> rankList;
	int index = 0;
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
			BuildMutationTemplate* mt = new BuildMutationTemplate(lines, bbLib, scLib, scLibSimp, atLib, rn);
			int designAA = mt->singleResidueDesign(ec, eS1S2);
			float sai = mt->targetNode->sai;
			int saiIndex = (int)(sai*4);
			int natAA = mt->targetNode->aaType;
			float deltaVol = rn->getVol(designAA) - rn->getVol(natAA);

			result[designAA] += 1.0;
			result[20] += 1.0;
			if(designAA == natAA){
				result[21] += 1.0;
			}

			result[22 + saiIndex] += deltaVol;
			sprintf(xx, "%d", index);
			string fileid = string(xx);

			index++;
			delete mt;
		}
	}
}



void getMeanRank(const string& datFile, ProEnergyTable* pe, ProS1S2Energy* eS1S2, DesignPara* dp, ofstream& out, const string& rotTag){

	ifstream file;
	file.open(datFile, ios::in);
	if(!file.is_open()){
		cout << "fail to open file: " << datFile << endl;
		exit(1);
	}

	ResName* rn = new ResName();
	cout << "init scLib" << endl;
	ResScRotamerLib* scLib = new ResScRotamerLib(rotTag);
	ResScRotamerLibMini* scLibSimp = new ResScRotamerLibMini();
	cout << "init bbLib" << endl;
	ResBBRotamerLib* bbLib = new ResBBRotamerLib();

	AtomLib* atLib = new AtomLib();
	cout << "init ec" << endl;
	ProParameter* pa;
	EnergyCalculator* ec = new EnergyCalculator(pe, dp, pa);

	cout << "finish ec" << endl;
	string s;
	vector<string> spt;
	vector<string> lines;
	//vector<int> rankList;
	int index = 0;
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

			BuildMutationTemplate* mt = new BuildMutationTemplate(lines, bbLib, scLib, scLibSimp, atLib, rn);
			//mt->printBBInfo();
			//string rankInfo = mt->nativeRank(ec, eS1S2);
			string type = rn->intToTri(mt->targetNode->aaType);

			//int rankS1 = mt->rankS1(eS1S2);
			//int rankS2 = mt->rankS2(eS1S2);
			//int rankS1S2 = mt->rankS1S2(eS1S2);
			//int rankAtomic = mt->rankAtomic(ec);
			string rankTotal = mt->nativeRank(ec, eS1S2);
			//int rankAtomicABACUS = mt->rankAtomicABACUS(ec);

			//int rankABACUS = mt->nativeRankABACUS(ec, eS1S2);
			//double logits = mt->nativeLogits(ec, eS1S2);
			//sprintf(xx, "%d", index);
			out << rankTotal << endl;
			//out << type << " " << rankS1 << " " << rankS2 << endl;

			//out << type << " " << rankS1 << " " << rankS2 << " " << rankS1S2 << " " << endl;
			//out << rankS1 << " " << rankS2 << " " << rankS1S2 << endl;
			//rankList.push_back(rank);
			index++;
			delete mt;
		}
	}


	delete rn;
	delete scLib;
	delete bbLib;
	delete ec;
}




void testUbq() {
	string datFile = "/user/xiongpeng/proteinModeling/resDesign/data/train/1ubq.dat";
	string outFile = "/user/xiongpeng/1.out";
	ProParameter* pa = new ProParameter();
	DesignPara* dp = new DesignPara("/user/xiongpeng/sword/designPara/p5");
	ofstream out;
	out.open(outFile.c_str(), ios::out);
	ProEnergyTable* pe = new ProEnergyTable(dp);
	ProS1S2Energy* eS1S2= new ProS1S2Energy(dp, 0);
	getMeanRank(datFile, pe, eS1S2, dp, out, "");
	out.close();
}


void checkBBInfo() {
	string datFile = "/user/xiongpeng/proteinModeling/resDesign/data/train/1ubq.dat";
	string outFile = "/user/xiongpeng/1.out";
	ProParameter* pa = new ProParameter();
	DesignPara* dp = new DesignPara("/user/xiongpeng/sword/designPara/p5");
	ofstream out;
	out.open(outFile.c_str(), ios::out);
	ProEnergyTable* pe = new ProEnergyTable(dp);
	ProS1S2Energy* eS1S2= new ProS1S2Energy(dp, 0);
	EnergyCalculator* ec = new EnergyCalculator(pe, dp, pa);


	string input = "/user/xiongpeng/input";
	//SeqDesignTemplate* dt = new SeqDesignTemplate(input, ec, eS1S2);
	//dt->printBackboneInfo();
	//delete dt;
}



void pdb400Test(const string& pdbID){
	string datFile = "/user/xiongpeng/sword/resDesign/" + pdbID + ".dat";
	string outFile = "/user/xiongpeng/sword/resDesign/sword/" + pdbID + ".out";
	ProParameter* pa = new ProParameter();
	DesignPara* dp = new DesignPara("/user/xiongpeng/sword/designPara/py");
	ofstream out;
	out.open(outFile.c_str(), ios::out);
	ProEnergyTable* pe = new ProEnergyTable(dp);
	ProS1S2Energy* eS1S2= new ProS1S2Energy(dp, 0);
	getMeanRank(datFile, pe, eS1S2, dp, out, "2");
	out.close();
}

void mutationFinal(const string& dat,  const string& designPara, const string& result, const string& rotTag){
	DesignPara* dp = new DesignPara("/user/xiongpeng/sword/designPara/" + designPara);
	ofstream out;
	out.open("/user/xiongpeng/proteinModeling/resDesign/result/" + result, ios::out);
	cout << "et" << endl;
	ProEnergyTable* pe = new ProEnergyTable(dp);
	cout << "s1s2" << endl;
	ProS1S2Energy* eS1S2= new ProS1S2Energy(dp, 0);
	cout << "get rank" << endl;
	getMeanRank(dat, pe, eS1S2, dp, out, rotTag);
	out.close();
}


void testRankFromDesignTemplate(const string& designPara, string& result) {
	ifstream file;
	file.open("/user/xiongpeng/sword/train/list3", ios::in);
	vector<string> pdbList;
	string s;
	while(getline(file,s)){
		pdbList.push_back(s);
	}

	cout << "init pa" << endl;
	ProParameter* pa = new ProParameter();
	DesignPara* dp = new DesignPara("/user/xiongpeng/sword/designPara/"+designPara);
	ofstream out;
	out.open("/user/xiongpeng/proteinModeling/resDesign/result/" + result, ios::out);

	cout << "init et" << endl;
	ProEnergyTable* pe = new ProEnergyTable();
	cout << "init s1s2" << endl;
	ProS1S2Energy* eS1S2= new ProS1S2Energy(dp, 0);
	cout << "init ec" << endl;
	EnergyCalculator* ec = new EnergyCalculator(pe, dp, pa);
	ResBBRotamerLib* bbLib = new ResBBRotamerLib();
	ResScRotamerLib* scLib = new ResScRotamerLib();
	AtomLib* atLib = new AtomLib();

	for(int i=0;i<pdbList.size();i++){
		string pdbID = pdbList[i];
		string inputFile = "/user/xiongpeng/sword/train/"+pdbID+"/input";
		string s2File = "/user/xiongpeng/sword/train/S2/"+pdbID+".s2";

		SeqDesignTemplate* dt = new SeqDesignTemplate(inputFile, bbLib, scLib, atLib, ec, eS1S2);
		dt->loadS1S2(s2File);
		dt->getPositionRanks(out);

		delete dt;
	}
	out.close();
}




int main(int argc, char** argv){

	//string dat = string(argv[1]);
	string para = string(argv[1]);
	string result = string(argv[2]);
	//string pdbID = string(argv[3]);

	testRankFromDesignTemplate(para, result);
	//checkBBInfo();

	//testUbq();


	//pdb400Test(string(argv[1]));

}


