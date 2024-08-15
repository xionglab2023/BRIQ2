#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "model/StructureModel.h"
#include "predNA/BackboneModelingTemplate.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredna;
using namespace std;


int main(int argc, char** argv){

	string pdbFile = string(argv[1]);
	string output = string(argv[2]);

    string cnt = "";
    vector<string> lines;

    ifstream file;
    file.open(pdbFile, ios::in);
    string line;
    while(getline(file, line)){
        lines.push_back(line);
        if(line.length() > 3 && line.substr(0,3)=="cnt")
            cnt = line.substr(5, line.length()-5);
    }


	
    RNAPDB* pdb = new RNAPDB(pdbFile);
	vector<RNABase*> baseList = pdb->getBaseList();
	

    int len = baseList.size();

    if(cnt.length() != len){
        cout << "base num: " << len << endl;
        cout << "cnt length: " << cnt.length() << endl;
        exit(0);
    }

    bool* tagList = new bool[len];
    for(int i=0;i<len;i++){
        tagList[i] = false;
    }

    int selectBaseNum = 0;
    cout << "cnt: " << cnt << endl;
    for(int i=0;i<len-1;i++){
        if(baseList[i]->connectToNeighbor(baseList[i+1]))
            cout << "yes" << endl;
        else 
            cout << "no" << endl;
        if(cnt[i] == '-' && !baseList[i]->connectToNeighbor(baseList[i+1])){
            selectBaseNum++;
            for(int j=i-2;j<=i+2;j++){
                if(j>=0 && j < len){
                    tagList[j] = true;
                }
            }
        }
    }

    

    if(selectBaseNum == 0){
	    delete pdb;
        ofstream out;
        out.open(output, ios::out);
        for(int i=0;i<lines.size();i++){
            out << lines[i] << endl;
        }
        out.close();
        exit(0);
    }
    else {
        char* csnList = new char[len];
        for(int i=0;i<len;i++){
            if(tagList[i])
                csnList[i] = '0';
            else 
                csnList[i] = 'F';
        }
        string csn = string(csnList);
        cout << "cnt: " << cnt << endl;
        cout << "csn: " << csn << endl;
	    ForceFieldPara* para = new ForceFieldPara();
	    BackboneModelingTemplate* bm = new BackboneModelingTemplate(pdb, cnt, csn, para);
	    double rms = bm->runMC();
	    BRTreeInfo* bi = bm->toTreeInfo();
	    bi->printPDB(output);
        delete [] csnList;
        delete pdb;
	    delete para;
	    delete bm;
    }

 




}