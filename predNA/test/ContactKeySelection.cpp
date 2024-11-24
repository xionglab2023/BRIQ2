#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include "model/StructureModel.h"
#include "predNA/NuGraph.h"
#include "predNA/NuSampling.h"
#include "tools/CmdArgs.h"
#include "tools/ThreadPool.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredNA;
using namespace std;
using namespace NSPtools;
using namespace NSPthread;


void printHelp(){
	cout << "Usage: " << endl;
	cout << "briqx_keySelection -in $INPUTFILE -key $KEYFILE -n $MAXNUM -out $OUTPUT" <<endl;
}

class keyEnergy{
    public:

    string key;
    double ene;

    keyEnergy(){
        this->key = "";
        this->ene = 0.0;
    }
    
    keyEnergy(const string& key, double ene){
        this->key = key;
        this->ene = ene;
    }
};


bool cmp_energy(keyEnergy* k1, keyEnergy* k2){
	return (k1->ene < k2->ene);
}

double keySimilarity(int seqLen, int* typeList, int* sepTable, const string& keyA, const string& keyB, BasePairLib* bpLib){

    if(keyA.length() % 6 != 0 || keyB.length() %6 != 0) {
        cout << "invalid key length" << endl;
        exit(0);
    }

    int bpMtx1[seqLen*seqLen];
    int bpMtx2[seqLen*seqLen];

    for(int i=0;i<seqLen*seqLen;i++){
        bpMtx1[i] = -1;
        bpMtx2[i] = -1;
    }

    for(int i=0;i<keyA.length();i+=6){
        int id1 = (keyA[i]-'!')*90 + (keyA[i+1] - '!');
        int id2 = (keyA[i+2]-'!')*90 + (keyA[i+3] - '!');
        int clusterID = (keyA[i+4]-'!')*90 + (keyA[i+5] - '!');
        bpMtx1[id1*seqLen+id2] = clusterID;
    }

    for(int i=0;i<keyB.length();i+=6){
        int id1 = (keyB[i]-'!')*90 + (keyB[i+1] - '!');
        int id2 = (keyB[i+2]-'!')*90 + (keyB[i+3] - '!');
        int clusterID = (keyB[i+4]-'!')*90 + (keyB[i+5] - '!');
        bpMtx2[id1*seqLen+id2] = clusterID;
      
    }

    double totalScore = 0.0;
    double simScore = 0.0;
    int clusterID1, clusterID2, sep;
    double ene1, ene2, minEne;

    for(int i=0;i<seqLen;i++){
        for(int j=i+1;j<seqLen;j++){
            clusterID1 = bpMtx1[i*seqLen+j];
            clusterID2 = bpMtx2[i*seqLen+j];
            sep = sepTable[i*seqLen+j];
            if(clusterID1 < 0 && clusterID2 < 0) continue;
            ene1 = 0.0;
            ene2 = 0.0;
            minEne = 0.0;


            if(clusterID1 >= 0)
                ene1 = bpLib->getEnergy(clusterID1, typeList[i], typeList[j], sep);
            if(clusterID2 >= 0)
                ene2 = bpLib->getEnergy(clusterID2, typeList[i], typeList[j], sep);

            minEne = ene1;
         
            if(ene2 < minEne) {
                minEne = ene2;
            }

            if(clusterID1 == clusterID2) {
                totalScore += minEne;
                simScore += minEne;
            }
            else if(clusterID1 >= 0 && clusterID2 >= 0){
                
                double d = bpLib->getDMDistance(clusterID1, clusterID2, typeList[i], typeList[j], sep);
                if(d > 2.0) {
                    totalScore += minEne;
                }
                else if(d < 1.5) {
                    totalScore += minEne;
                    simScore += minEne;
                }
                else {
                    totalScore += minEne;
                    simScore += minEne * (2.0 - d) * 2.0;
                }
            }
            else {
                totalScore += minEne;
            }
        }
    }

    double p =  1.0*simScore/totalScore;  

    return p;  
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

    string inputFile = cmdArgs.getValue("-in");
    string keyFile = cmdArgs.getValue("-key");
    string outputFile = cmdArgs.getValue("-out");

    int maxNum = 100;

    if(cmdArgs.specifiedOption("-n"))
        maxNum = atoi(cmdArgs.getValue("-n").c_str());

    ifstream file;
    file.open(keyFile.c_str(), ios::in);
    if(!file.is_open()){
        cout << "fail to open file: " << inputFile << endl;
        exit(0);
    }
    
    string key;
    double ene;
    vector<string> keys;
    vector<double> eneList;

    string s;
    vector<string> spt;

    while(getline(file, s)){
        splitString(s, " ", &spt);
        if(spt.size() < 2) continue;
        keys.push_back(spt[0]);
        eneList.push_back(atof(spt[1].c_str()));    
    }

    file.close();

    NSPtools::InputParser input(inputFile);
    string seq = input.getValue("seq");
	string cnt = input.getValue("cnt");
    int seqLen = seq.length();

    cout << seq << endl;
    cout << cnt << endl;

    int typeList[seqLen];
    bool connectToDownstream[seqLen];
	int sepTable[seqLen*seqLen];

    int i,j,k;
	/*
	 * read chain break points
	 */

    map<char, int> typeMap;
    typeMap['A'] = 0;
    typeMap['U'] = 1;
    typeMap['G'] = 2;
    typeMap['C'] = 3;
    typeMap['a'] = 4;
    typeMap['t'] = 5;
    typeMap['g'] = 6;
    typeMap['c'] = 7;

    for(i=0;i<seqLen;i++){
        typeList[i] = typeMap[seq[i]];
    }

	for(i=0;i<seqLen;i++){
		connectToDownstream[i] = true;
	}
	connectToDownstream[seqLen-1] = false;

	for(i=0;i<seqLen-1;i++){
		if(cnt[i] == '|'){
			connectToDownstream[i] = false;
		}
		else if(cnt[i] == '-'){
			connectToDownstream[i] = true;
		}
		else {
			cout << "invalid cnt string: " << cnt << endl;
			cout << "cnt example: ------|-------" << endl;
			exit(0);
		}
	}

	for(i=0;i<seqLen;i++){
		for(j=0;j<seqLen;j++){
			int ij = i*seqLen+j;
			if(i==j) sepTable[ij] = 0;
			else if(j == i+1 && connectToDownstream[i]) sepTable[ij] = 1;
			else if(j == i-1 && connectToDownstream[j]) sepTable[ij] = -1;
			else sepTable[ij] = 2;
		}
	}

    int n = keys.size();

    keyEnergy* keList[n];
    for(i=0;i<n;i++){
        keList[i] = new keyEnergy(keys[i], eneList[i]);
    }



    BasePairLib* bpLib = new BasePairLib("stat");

    

    vector<string> selectKeys;
    vector<double> selectEnes;



    sort(keList, keList+n, cmp_energy);
    bool selected[n];
    



    double distCutoff1 = 1.0;
    double distCutoff2 = 1.0;
    double distCutoff3 = 1.0;

    double eneCutoff1 = 0.8*keList[0]->ene;
    double eneCutoff2 = 0.6*keList[0]->ene;
    double eneCutoff3 = 0.4*keList[0]->ene;


    if(keList[0]->ene > 0) {
        eneCutoff1 = 1.5*keList[0]->ene;
        eneCutoff2 = 2.0*keList[0]->ene;
        eneCutoff3 = 3.0*keList[0]->ene;
    }

    double distCutoff;

    vector<double> cutoffList;


    double stepLength = -0.02;
    int round = 0;

    while(selectKeys.size()==0 ||  selectKeys.size() > maxNum*1.05 || selectKeys.size() < maxNum ){
        round++;
        if(round > 20) break;

        if(selectKeys.size() > maxNum*1.1 && stepLength >0){
            stepLength = -stepLength*0.7;
        }
        else if(selectKeys.size() > 0 && selectKeys.size() < maxNum && stepLength <0){
            stepLength = -stepLength*0.7;
        }

        if(abs(stepLength) < 0.001) break;

        distCutoff1 = distCutoff1 + stepLength;
        distCutoff2 = 1 - (1-distCutoff1)*2;
        distCutoff3 = 1 - (1-distCutoff1)*4;

        selectKeys.clear();
        selectEnes.clear();
        for(i=0;i<n;i++){
            selected[i] = false;
        }

        for(i=0;i<n;i++){
            if(selected[i]) continue;
            double ene = keList[i]->ene;

            if(ene > eneCutoff3) continue;
            else if(ene > eneCutoff2) {
                distCutoff = distCutoff3;
            }
            else if(ene > eneCutoff1) {
                distCutoff = distCutoff2;
            }
            else 
                distCutoff = distCutoff1;


            selectKeys.push_back(keList[i]->key);
            selectEnes.push_back(keList[i]->ene);

            for(j=i+1;j<n;j++){
                double d =  keySimilarity(seqLen, typeList, sepTable, keList[i]->key, keList[j]->key, bpLib);
                if(d > distCutoff){
                    selected[j] = true;
                }
            }
        }

    }

    int maxLen = 200;
    for(i=0;i<selectKeys.size();i++){
        if(selectKeys[i].length() > maxLen){
            maxLen = selectKeys[i].length();
        }
    }

    ofstream out;
    out.open(outputFile.c_str(), ios::out);
    char xx[maxLen + 100];


    for(i=0;i<selectKeys.size() && i < maxNum;i++){
    
        sprintf(xx, "%s %7.3f", selectKeys[i].c_str(), selectEnes[i]);
        out << string(xx) << endl;
    }
   
    out.close();

    delete bpLib;

}