#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include "model/StructureModel.h"
#include "predNA/NuGraph.h"
#include "tools/CmdArgs.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredNA;
using namespace std;
using namespace NSPtools;

class pdbEne{
	public:
	RNAPDB* pdb;
	double ene;

	pdbEne(RNAPDB* pdb, double ene){
		this->pdb = pdb;
		this->ene = ene;
	}
};

bool pdbEneCmp(pdbEne* peA, pdbEne* peB){
	return peA->ene < peB->ene;
}

void selectModel(vector<RNAPDB*>& pdbList, int n, double rmsdCutoff, const string& output){

	int totalNum = pdbList.size();
	cout << "totalNum: " << totalNum << endl;

	if(totalNum == 0) 
		return;

	pdbEne* peList[totalNum];
	for(int i=0;i<pdbList.size();i++){
		peList[i] = new pdbEne(pdbList[i], pdbList[i]->ene);
	}
	
    sort(peList, peList+totalNum, pdbEneCmp);

	double minEne = peList[0]->ene;
	double eneCutoff = minEne * 0.9;
	char xx[200];
	ofstream out;


	if(minEne >= 0) {
		cout << "positive energy, select only one model: " << endl;
		for(int i=0;i<1;i++){
			if(output.length() > 3 && output.substr(output.length()-3, 3) == "pdb")
				sprintf(xx, "%s-%d.pdb", output.substr(0, output.length()-3), i);
			else 
				sprintf(xx, "%s-%d.pdb", output, i);
		
			string outputFile = string(xx);
			out.open(outputFile.c_str());
			peList[i]->pdb->printPDBFormat(out);
			out.close();
		}
	}

	vector<pdbEne*> selectList;
	

	int round = 0;
	int i,j;
	while(true){

		selectList.clear();
		bool taged[totalNum];
		for(i=0;i<totalNum;i++){
			if(peList[i]->ene > eneCutoff) //ignore higher energy model
				taged[i] = true;

            taged[i] = false;
        }

        for(i=0;i<totalNum;i++){
			

            if(taged[i]) continue;
			selectList.push_back(peList[i]);
            for(j=i+1;j<totalNum;j++){
				if(taged[j]) continue;
                double d = peList[i]->pdb->baseRMSD(peList[j]->pdb);
                if(d < rmsdCutoff){
                    taged[j] = true;
                }
            }
        }

		round++;

		int selectNum = selectList.size();
		if(selectNum >= n || round > 20) break;
		eneCutoff = eneCutoff - minEne * 0.05;
	}


	for(int i=0;i<n && i<selectList.size();i++){
		if(output.length() > 3 && output.substr(output.length()-4, 4) == ".pdb")
			sprintf(xx, "%s-%d.pdb", output.substr(0, output.length()-4).c_str(), i);
		else 
			sprintf(xx, "%s-%d.pdb", output.c_str(), i);
		
		string outputFile = string(xx);
		out.open(outputFile.c_str());
		selectList[i]->pdb->printPDBFormat(out);
		out << "ene " << selectList[i]->ene << endl;
		out.close();
	}

	for(int i=0;i<n && i<selectList.size();i++){
		RNAPDB* pdb1 = selectList[i]->pdb;
		for(int j=0;j<n && j < selectList.size();j++){
			RNAPDB* pdb2 = selectList[j]->pdb;
			double rms = pdb1->baseRMSD(pdb2);
			printf("%6.3f ",  rms);
		}
		cout << endl;
	}
}

void readPDBFromTreeFile(const string& treeFile, vector<RNAPDB*>& pdbList){

	int seqLen = 0;

	ifstream file;
	file.open(treeFile.c_str(), ios::in);

	string s;
	vector<string> spt;
	string seq = "";
	string cnt = "";
	double ene;
	vector<int> rotIDList;
	vector<double> d1List;
	vector<double> d2List;
	vector<CsMove> cmList;
	vector<NuNode*> nodes;

	map<char,int> augcMap;
	augcMap['A'] = 0;
	augcMap['U'] = 1;
	augcMap['G'] = 2;
	augcMap['C'] = 3;
	augcMap['a'] = 4;
	augcMap['t'] = 5;
	augcMap['g'] = 6;
	augcMap['c'] = 7;

	AtomLib* atLib = new AtomLib();
	ForceFieldPara* para = new ForceFieldPara();
	BaseRotamerLib* baseRotLib = new BaseRotamerLib(atLib);
	RiboseRotamerLib* riboLib = new RiboseRotamerLib(para);

	int i,j;
	LocalFrame cs0;
	LocalFrame cs1;

	
	string augc = "AUGCatgc";
	char ss[20];
	int seqID = 0;

	while(getline(file, s)){
		splitString(s, " ", &spt);
		if(spt[0] == "model"){
			rotIDList.clear();
			d1List.clear();
			d2List.clear();
			cmList.clear();
			nodes.clear();

			seq = spt[1];
			cnt = spt[2];
			ene = atof(spt[3].c_str());
			seqLen = seq.length();
			getline(file, s);
			splitString(s, " ", &spt);
			for(i=0;i<seqLen;i++){
				rotIDList.push_back(atoi(spt[i].c_str()));
			}
			getline(file, s);
			splitString(s, " ", &spt);
			for(i=0;i<seqLen;i++){
				d1List.push_back(atof(spt[i].c_str()));
			}
			getline(file, s);
			splitString(s, " ", &spt);
			for(i=0;i<seqLen;i++){
				d2List.push_back(atof(spt[i].c_str()));
			}

			for(i=0;i<seqLen-1;i++){
				getline(file, s);
				CsMove cm(s);
				cmList.push_back(cm);
			}

			for(i=0;i<seqLen;i++){

				cs1 = LocalFrame();
				if(i > 0)
					cs1 = cs0 + cmList[i-1];

				int baseType = augcMap[seq[i]];
				int rotID = rotIDList[i];
				NuNode* node = new NuNode(i, baseType, cs1, baseRotLib->baseLib[baseType], riboLib->rotLib[baseType][rotID], atLib);
				
				if(cnt[i] == '-') {
					PhosphateRotamer* prot = new PhosphateRotamer(d1List[i], d2List[i]);
					node->phoConf->updateLocalFrameAndRotamer(node->riboseConf->cs2, prot, 0.0);	
				}
				
				nodes.push_back(node);
			}

			RNAChain* rc = new RNAChain("A");
			seqID = 0;
			for(i=0;i<seqLen;i++) {
				seqID++;
				sprintf(ss, "%d", seqID);
				RNABase* base = new RNABase(string(ss), "A", seq[i]);
				vector<Atom*> aList = nodes[i]->toAtomList(atLib);
				for(int j=0;j<aList.size();j++)
					base->addAtom(aList[j]);
				rc->addBase(base);
			}


			for(i=0;i<seqLen-1;i++){
				if(cnt[i] == '|') continue;
				vector<Atom*> aList = nodes[i]->toPhoAtomList(atLib);
				RNABase* base = rc->getBaseList()[i+1];
				for(j=0;j<aList.size();j++)
					base->addAtom(aList[j]);
			}

			RNAPDB* pdb = new RNAPDB();
			pdb->addChain(rc);
			pdb->ene = ene;
			pdbList.push_back(pdb);

		}
	}
}

void printHelp(){
	cout << "Usage: " << endl;
	cout << "briqx_modelSelection -pdbList $PDBLIST -n $MODELNUM -cutoff $RMSDCUTOFF -out $OUTPUT" <<endl;
	cout << "briqx_modelSelection -tre $TREEFILE -n $MODELNUM -cutoff $RMSDCUTOFF -out $OUTPUT" << endl;
	cout << "briqx_modelSelection -treList $TREELIST -n $MODELNUM -cutoff $RMSDCUTOFF -out $OUTPUT" << endl;
}

int main(int argc, char** argv){

    
    if(argc == 1) {
        printHelp();
        return EXIT_SUCCESS;
    }

	clock_t start = clock();
    CmdArgs cmdArgs{argc, argv};
	if(cmdArgs.specifiedOption("-h")){
		printHelp();
        return EXIT_SUCCESS;
	}

	int n = 10;

	if(cmdArgs.specifiedOption("-n"))
    	n = atoi(cmdArgs.getValue("-n").c_str());

	double cutoff = 2.0;

	if(cmdArgs.specifiedOption("-cutoff"))
		cutoff = atof(cmdArgs.getValue("-cutoff").c_str());


	string output = cmdArgs.getValue("-out");

	vector<RNAPDB*> pdbList;
	vector<double> eneList;

	string s;

	if(cmdArgs.specifiedOption("-plot")) {
		string tempFile = cmdArgs.getValue("-temp");
		RNAPDB* temp = new RNAPDB(tempFile);

		if(cmdArgs.specifiedOption("-tre")){
			readPDBFromTreeFile(cmdArgs.getValue("-tre"), pdbList);
		}
		else if(cmdArgs.specifiedOption("-treList")){
			ifstream file;
			file.open(cmdArgs.getValue("-treList").c_str(), ios::in);
			if(!file.is_open()){
				cout << "fail to open file: " << cmdArgs.getValue("-treList") << endl;
			}
			while(getline(file, s)){
				readPDBFromTreeFile(s, pdbList);
			}
			file.close();
		}
		else {
			cout << "please specify tre file" << endl;
			exit(1);
		}

		ofstream f;
		f.open(output, ios::out);

		for(int i=0;i<pdbList.size();i++){
			double ene = pdbList[i]->ene;
			double rms = pdbList[i]->baseRMSD(temp);
			f << rms << " " << ene << endl;
		}
		f.close();

		delete temp;

	}
	else if(cmdArgs.specifiedOption("-pdbList")) {
		ifstream file;
		file.open(cmdArgs.getValue("-pdbList").c_str(), ios::in);
		if(!file.is_open()){
			cout << "fail to open file: " << cmdArgs.getValue("-pdbList") << endl;
		}
		while(getline(file, s)){
			RNAPDB* pdb = new RNAPDB(s);
			pdbList.push_back(pdb);
		}
		file.close();
		selectModel(pdbList, n, cutoff, output);
	}
	else if(cmdArgs.specifiedOption("-tre")) {
		
		readPDBFromTreeFile(cmdArgs.getValue("-tre"), pdbList);
		selectModel(pdbList, n, cutoff, output);
		
	}
	else if(cmdArgs.specifiedOption("-treList")) {
		ifstream file;
		file.open(cmdArgs.getValue("-treList").c_str(), ios::in);
		if(!file.is_open()){
			cout << "fail to open file: " << cmdArgs.getValue("-treList") << endl;
		}
		while(getline(file, s)){
			readPDBFromTreeFile(s, pdbList);
		}
		file.close();
		selectModel(pdbList, n, cutoff, output);
	}
	else {
		printHelp();
		exit(1);
	}



}
