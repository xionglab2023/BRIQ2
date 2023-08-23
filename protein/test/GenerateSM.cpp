/*
 * GenerateSM.cpp
 *
 */

#include "model/StructureModel.h"
#include "model/ResBBRotamerLib.h"
#include "math/AAProbabilityMatrix.h"
#include "math/AAScoreMatrix.h"
#include "math/AAProbabilityArray.h"
#include "math/AAScoreArray.h"
#include "math/Stat.h"
#include "protein/SeqDesignTemplate.h"
#include "para/ProParameter.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>

using namespace std;
using namespace NSPmodel;
using namespace NSPmath;
using namespace NSPforcefield;
using namespace NSPprotein;

void generateSM(ResPairInfo* pi, ResBBRotamerLib* bbLib, ofstream& out, double bwSai, double initWt, double sumCutoff){
	vector<string> fileList;
	char xx[200];
	for(int i=0;i<10;i++){
		for(int j=0;j<10;j++){
			double disA = i*0.1+0.05 - pi->saiA;
			double disB = j*0.1+0.05 - pi->saiB;
			double d = sqrt(disA*disA + disB*disB);
			if(d > 0.45) continue;
			int k = i*10 + j;
			sprintf(xx, "/user/xiongpeng/sword/data50/split/saiSub/%s-%d", pi->key.c_str(), k);
			fileList.push_back(string(xx));
		}
	}

	ResName rn;

	int fileNum = fileList.size();
	int n = 0;

	ifstream f;
	string line;

	vector<float> d1List;
	vector<float> d2List;
	vector<int> aaList1;
	vector<int> aaList2;

	map<string, double> bwMap;
	bwMap["CC1"] =  0.0863;
	bwMap["CC2"] =  0.3050;
	bwMap["CC3"] =  0.4180;
	bwMap["CC4"] =  0.4590;
	bwMap["CC5"] =  0.7160;
	bwMap["CE1"] =  0.1937;
	bwMap["CE2"] =  0.2820;
	bwMap["CE3"] =  0.3830;
	bwMap["CE4"] =  0.4140;
	bwMap["CE5"] =  0.5960;
	bwMap["CH1"] =  0.1755;
	bwMap["CH2"] =  0.2340;
	bwMap["CH3"] =  0.3020;
	bwMap["CH4"] =  0.3190;
	bwMap["CH5"] =  0.6700;
	bwMap["EC1"] =  0.1863;
	bwMap["EC2"] =  0.2870;
	bwMap["EC3"] =  0.3930;
	bwMap["EC4"] =  0.3680;
	bwMap["EC5"] =  0.5710;
	bwMap["EE1"] =  0.0899;
	bwMap["EE2"] =  0.2330;
	bwMap["EE3"] =  0.3410;
	bwMap["EE4"] =  0.3350;
	bwMap["EE5"] =  0.3320;
	bwMap["EH1"] =  0.1501;
	bwMap["EH2"] =  0.2720;
	bwMap["EH3"] =  0.3160;
	bwMap["EH4"] =  0.3830;
	bwMap["EH5"] =  0.5560;
	bwMap["HC1"] =  0.1164;
	bwMap["HC2"] =  0.1380;
	bwMap["HC3"] =  0.1850;
	bwMap["HC4"] =  0.1920;
	bwMap["HC5"] =  0.6670;
	bwMap["HE1"] =  0.0977;
	bwMap["HE2"] =  0.2010;
	bwMap["HE3"] =  0.2390;
	bwMap["HE4"] =  0.2940;
	bwMap["HE5"] =  0.5590;
	bwMap["HH1"] =  0.0349;
	bwMap["HH2"] =  0.0669;
	bwMap["HH3"] =  0.1009;
	bwMap["HH4"] =  0.0993;
	bwMap["HH5"] =  0.6000;

	double tmpBW1 = 1/bwSai;
	double tmpBW2 = 1/bwMap[pi->key]/initWt;

	ResPairInfo* q;
	double saiA, saiB, d1, d2;
	int bbA, bbB;
	string triA, triB;
	CsMove cm;
	DistanceMatrixRes dm;
	for(int i=0;i<fileNum;i++){
		f.open(fileList[i].c_str(), ios::in);
		if(!f.is_open()){
			cout << "fail to open " << fileList[i] << endl;
			continue;
		}
		while(getline(f, line)){
			saiA = atof(line.substr(6,5).c_str());
			saiB = atof(line.substr(23, 5).c_str());
			d1 = pi->saiA - saiA;
			d2 = pi->saiB - saiB;
			triA = line.substr(0,3);
			triB = line.substr(17,3);
			cm = CsMove(line.substr(34, 130));
			dm = DistanceMatrixRes(cm);

			d1List.push_back((float)sqrt(d1*d1+d2*d2));
			d2List.push_back((float)pi->dm.distanceTo(dm));
			aaList1.push_back(rn.triToInt(triA));
			aaList2.push_back(rn.triToInt(triB));
		}
		f.close();
	}

	double u1, u2, u3, u4;
	double sum = 0;
	double pa1[20];
	double pa2[20];
	double pm[400];
	double pmBg[400];

	while(sum < sumCutoff && tmpBW2 > 0.5){
		sum = 0;
		for(int i=0;i<20;i++) {
			pa1[i] = 0;
			pa2[i] = 0;
		}
		for(int i=0;i<20;i++) {
			for(int j=0;j<20;j++) {
				pm[i*20+j] = 0;
				pmBg[i*20+j] = 0;
			}
		}

		for(int k=0;k<d1List.size();k++) {
			u1 = d1List[k]*tmpBW1;
			u2 = d2List[k]*tmpBW2;
			double wt = exp(-0.5*(u1*u1+u2*u2));

			pa1[aaList1[k]] += wt;
			pa2[aaList2[k]] += wt;
			pm[aaList1[k]*20+aaList2[k]] += wt;
		}


		for(int j=0;j<20;j++) {
			for(int k=0;k<20;k++) {
				sum += pm[j*20+k];
			}
		}
		for(int i=0;i<20;i++) {
			for(int j=0;j<20;j++) {
				pmBg[i*20+j] = pa1[i]*pa2[j];
			}
		}

		tmpBW2 = tmpBW2/1.1;
		printf("%5.3f %8.2f\n", 1/tmpBW2, sum);
	}

	AAProbabilityMatrix* pm1 = new AAProbabilityMatrix(pm, sum);
	AAProbabilityMatrix* pm2 = new AAProbabilityMatrix(pmBg, sum);
	pm1->addPseudoCount(pm2, 100);
	AAScoreMatrix* sm = new AAScoreMatrix(pm1, pm2);
	sprintf(xx, "%d %d\n", pi->indexA, pi->indexB);
	out << string(xx);
	sm->print(out);

	delete pm1;
	delete pm2;
	delete sm;
}

void getSM(const string& pdbID) {
	ProParameter* pa = new ProParameter();
	DesignPara* dp = new DesignPara("/user/xiongpeng/sword/designPara/py");
	ProEnergyTable* pe = new ProEnergyTable();
	ProS1S2Energy* eS1S2= new ProS1S2Energy(dp, 0);

	EnergyCalculator* ec = new EnergyCalculator(pe, dp, pa);
	string inputFile = "/user/xiongpeng/sword/train/" + pdbID + "/input";
	string outfile = "/user/xiongpeng/sword/train/S2/"+pdbID+".s2";
	ResBBRotamerLib* bbLib = new ResBBRotamerLib();
	ResScRotamerLib* scLib = new ResScRotamerLib();
	AtomLib* atLib = new AtomLib();
	SeqDesignTemplate* dt = new SeqDesignTemplate(inputFile, bbLib, scLib, atLib, ec, eS1S2);
	dt->printS2(outfile);
}


int main(int argc, char** argv){

	getSM(string(argv[1]));

	/*
	string pdbID = string(argv[1]);
	string paraIndex = string(argv[2]);
	double bwSai = atof(argv[3]);
	double topbb = atof(argv[4]);
	double sumCutoff = atof(argv[5]);

	cout << "init bbLib" << endl;

	ResBBRotamerLib* bbLib = new ResBBRotamerLib();

	cout << "fi" << endl;

	string outFile = "/user/xiongpeng/sword/train/S2/" + pdbID + "-" + paraIndex + ".s2";
	string pairFile = "/user/xiongpeng/sword/train/pairs/" + pdbID +  ".rp";
	ifstream f;
	ofstream out;
	out.open(outFile.c_str(), ios::out);
	f.open(pairFile.c_str(), ios::in);
	if(!f.is_open()){
		cout << "fail to open file: " << pairFile << endl;
	}
	vector<string> lines;
	string line;
	while((getline(f, line))){
		lines.push_back(line);
	}
	vector<string> spt;
	for(int i=0;i<lines.size();i+=2){
		string s1 = lines[i];
		string s2 = lines[i+1];
		cout << s1 << endl;
		ResPairInfo rp(s2);
		splitString(s1, " ", &spt);
		int id1 = atoi(spt[0].c_str());
		int id2 = atoi(spt[1].c_str());
		rp.indexA = id1;
		rp.indexB = id2;
		generateSM(&rp, bbLib, out, bwSai, topbb, sumCutoff);
	}
	*/
}



