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

string pdbToKey(RNAPDB* pdb, BasePairLib* bpLib){

    vector<RNABase*> baseList = pdb->getBaseList();
	string key;
	int i,j, k, sep, clusterID;
	char a,b;
	double ene;

    int seqLen = baseList.size();
    int sepTable[seqLen*seqLen];
    bool connectToDownstream[seqLen];

    for(i=0;i<seqLen;i++){
		connectToDownstream[i] = true;
	}
	connectToDownstream[seqLen-1] = false;



		for(i=0;i<seqLen-1;i++){
			if(baseList[i]->connectToNeighbor(baseList[i+1]))
				connectToDownstream[i] = true;
			else {
				connectToDownstream[i] = false;
			}
		}
	

    for(int i=0;i<seqLen;i++){
        for(int j=0;j<seqLen;j++){
            sepTable[i*seqLen+j] = 0;
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

	for(i=0;i<seqLen;i++){
		for(j=i+1;j<seqLen;j++){
			sep = sepTable[i*seqLen+j];
			double wt = 0.0;

			{
				BaseDistanceMatrix dm(*baseList[i], *baseList[j]);
				clusterID = bpLib->getPairType(dm, baseList[i]->baseTypeInt, baseList[j]->baseTypeInt, sep);
				if(clusterID < 0) continue;
                ene = bpLib->getEnergy(clusterID, baseList[i]->baseTypeInt, baseList[j]->baseTypeInt, sep);
                if(ene > -2.0) continue;

				a = i/90 + '!';
				b = i%90 + '!';
				key.push_back(a);
				key.push_back(b);

				a = j/90 + '!';
				b = j%90 + '!';
				key.push_back(a);
				key.push_back(b);

				a = clusterID/90 + '!';
				b = clusterID%90 + '!';
				key.push_back(a);
				key.push_back(b);
			}		
		}
	}
	return key;

}


double keyCoverage(int seqLen, int* typeList, int* sepTable, vector<string>& keyList, const string& natKey, BasePairLib* bpLib, bool printBpMtx){


    int bpMtxNat[seqLen*seqLen];
    int count[seqLen*seqLen];

    for(int i=0;i<seqLen*seqLen;i++){
        bpMtxNat[i] = -1;
        count[i] = 0;
    }

    for(int i=0;i<natKey.length();i+=6){
        int id1 = (natKey[i]-'!')*90 + (natKey[i+1] - '!');
        int id2 = (natKey[i+2]-'!')*90 + (natKey[i+3] - '!');
        int clusterID = (natKey[i+4]-'!')*90 + (natKey[i+5] - '!');
        bpMtxNat[id1*seqLen+id2] = clusterID;
    }

    for(int i=0;i<keyList.size();i++){
        for(int j=0;j<keyList[i].length();j+=6){
            int id1 = (keyList[i][j]-'!')*90 + (keyList[i][j+1] - '!');
            int id2 = (keyList[i][j+2]-'!')*90 + (keyList[i][j+3] - '!');
            int clusterID = (keyList[i][j+4]-'!')*90 + (keyList[i][j+5] - '!');
            if(bpMtxNat[id1*seqLen+id2] == clusterID) {
                count[id1*seqLen+id2] ++;
            }
        }
    }

    string augc = "AUGC";



    double totalScore = 0.0;
    double simScore = 0.0;
    int clusterID1, clusterID2, sep, n;

    double eneNat = 0.0, enePred = 0.0;

    for(int i=0;i<seqLen;i++){
        for(int j=i+1;j<seqLen;j++){
            clusterID1 = bpMtxNat[i*seqLen+j];
            n = count[i*seqLen+j];

            sep = sepTable[i*seqLen+j];
            if(clusterID1 < 0) continue;
            double ene = bpLib->getEnergy(clusterID1, typeList[i], typeList[j], sep);
            eneNat += ene;
            if(n>0) {
                enePred += ene;
            }

            if(printBpMtx) {
                char ta = augc[typeList[i]];
                char tb = augc[typeList[j]];
                printf("%-2d %-2d %c%c-%-4d %7.3f %4d\n",i, j, ta, tb, clusterID1, ene, n);
            }
        }
    }



    double p =  1.0*enePred/eneNat;  

    return p;  
}

void printHelp(){
	cout << "Usage: " << endl;
	cout << "briqx_keyAnalysis -in $INPUTFILE -key $KEYFILE -pdb $PDBFILE" <<endl;
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
    string natPDB = cmdArgs.getValue("-pdb");

    cout << "read pdb file" << endl;
    RNAPDB* pdb = new RNAPDB(natPDB);


    cout << "read key file: " << endl;
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

    cout << "read input file" << endl;

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

    cout << "init bpLib" << endl;
    BasePairLib* bpLib = new BasePairLib("stat");


    string natKey = pdbToKey(pdb, bpLib);
    cout << "nat Key: " << natKey << endl;



    vector<string> keyList;
    for(i=0;i<n;i++){
        vector<string> singleKey;
        singleKey.push_back(keys[i]);
        keyList.push_back(keys[i]);
        double coverage1 = keyCoverage(seqLen, typeList, sepTable, singleKey, natKey, bpLib,false);
        double coverage2 = keyCoverage(seqLen, typeList, sepTable, keyList, natKey, bpLib, false);

        printf("%-3d %6.4f %6.4f\n", i, coverage1, coverage2);
    }

    keyCoverage(seqLen, typeList, sepTable, keyList, natKey, bpLib, true);
}