

#include <time.h>
#include <iostream>
#include "forcefield/ProS1S2Energy.h"
#include "model/StructureModel.h"
#include "model/AssignProSSAndSasa.h"
#include "model/ResBBRotamerLib.h"
#include "model/ResBBRotamer.h"

using namespace std;
using namespace NSPforcefield;
using namespace NSPmodel;

void testS2(){
	DesignPara* para = new DesignPara("filename");
	ProS1S2Energy* et = new ProS1S2Energy(para, 0);
	string rpInfo = "PRO C 0.350  629 GLY H 0.350  683 -1.176543  1.381442   4.645712   -0.7438233 0.4958244  -0.4482021 -0.1152376 0.5654040  0.8167243  0.6583671  0.6591483  -0.3634231 5   1";
	ResPairInfo rp(rpInfo);

	AAScoreMatrix sm1 = et->getS2(&rp);
	AAScoreMatrix sm2 = et->getS2NearestNb(&rp);

	for(int i=0;i<20;i++) {
		for(int j=0;j<20;j++){
			printf("%-7.4f ", sm2.sm[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	for(int i=0;i<20;i++) {
		for(int j=0;j<20;j++){
			printf("%-7.4f ", sm1.sm[i][j]);
		}
		printf("\n");
	}
}

void testS1S2(const string& pdbFile){
	cout << "testS1S2: " << endl;
	clock_t start = clock();
	cout << "init S1S2 table: " << endl;

	string designPara = "/user/xiongpeng/sword/designPara/p5";

	DesignPara* para = new DesignPara(designPara);
	ProS1S2Energy* et = new ProS1S2Energy(para, 0);

	PDB* pdb = new PDB(pdbFile, "xxxx");
	AssignProSSAndSasa* ssa = new AssignProSSAndSasa(pdb);
	ResSasaPoints* rsp = new ResSasaPoints();
	ssa->updateSS();
	ssa->updateSasa(rsp);

	vector<Residue*> resList = pdb->getResList();

	cout << "res num: " << resList.size()<< endl;

	vector<ResInfo*> riList;
	vector<ResPairInfo*> rpList;
	ResBBRotamerLib* bbLib = new ResBBRotamerLib();
	AtomLib* atLib = new AtomLib();

	for(int i=0;i<resList.size();i++){
		int aa = resList[i]->intName;
		char ss = ssa->getSS(i);
		float sai = ssa->getSASA(i);
		ResBBRotamer* rot;

		if(i == 0)
			rot = new ResBBRotamer(resList[i], atLib);
		else
			rot = new ResBBRotamer(resList[i-1], resList[i], atLib);

		int bbIndex = bbLib->getRotamerIndex1K(rot);
		ResInfo* ri = new ResInfo(aa, ss, sai, bbIndex);
		riList.push_back(ri);
	}

	int n = resList.size();
	for(int i=0;i<n;i++){
		char ssA = riList[i]->ss;
		float saiA = riList[i]->sai;
		int bbA = riList[i]->bbIndex;
		LocalFrame csA = resList[i]->coordSys;
		XYZ tbA = resList[i]->getCbCoord();
		for(int j=i+1;j<n;j++){
			char ssB = riList[j]->ss;
			float saiB = riList[j]->sai;
			int bbB = riList[j]->bbIndex;
			LocalFrame csB = resList[j]->coordSys;
			XYZ tbB = resList[j]->getCbCoord();
			if(tbA.distance(tbB) > 8.0) continue;
			int sep = j-i;
			if(sep > 5) sep = 5;
			CsMove cm = csB - csA;
			ResPairInfo* rp = new ResPairInfo(i,j,resList[i]->intName, resList[j]->intName, ssA, ssB, saiA, saiB, cm ,sep);
			cout << rp->key << endl;

			rpList.push_back(rp);
		}
	}

	for(int i=0;i<riList.size();i++){
		AAScoreArray sa = et->getS1(riList[i]);
		printf("%-2d %6.3f\n", i, sa.sa[riList[i]->aaType]);
	}
	cout << endl;

	for(int i=0;i<rpList.size();i++){
		AAScoreMatrix sm = et->getS2(rpList[i]);
		int idA = rpList[i]->indexA;
		int idB = rpList[i]->indexB;
		int typeA = riList[idA]->aaType;
		int typeB = riList[idB]->aaType;
		printf("%-2d %-2d %7.3f\n", idA, idB, sm.sm[typeA][typeB]);
	}
}

int main(int argc, char** argv){
	testS1S2(string(argv[1]));
}


