#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include "model/StructureModel.h"
#include "predNA/NuGraph.h"
#include "tools/CmdArgs.h"
#include "tools/ThreadPool.h"
#include "predNA/NuSampling.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredNA;
using namespace std;
using namespace NSPtools;
using namespace NSPthread;

int main(int argc, char** argv){

    CmdArgs cmdArgs{argc, argv};

    string inputFile = cmdArgs.getValue("-in");
    string outputFile = cmdArgs.getValue("-out");
    string keyFile = cmdArgs.getValue("-key");

    string libType = "stat";
    if(cmdArgs.specifiedOption("-lib")) {
        libType = cmdArgs.getValue("-lib");
    }

    ifstream file;
    file.open(keyFile.c_str(), ios::in);
    vector<string> keys;
    vector<double> eneList;
    string line;
    vector<string> spt;
    while(getline(file, line)){
        splitString(line, " ", &spt);
        keys.push_back(spt[0]);
        eneList.push_back(atof(spt[1].c_str()));
    }

    vector<double> accuracyList;

	BasePairLib* pairLib = new BasePairLib(libType);
	RotamerLib* rotLib = new RotamerLib();
    AtomLib* atLib = new AtomLib();
    NuGraph* graph = new NuGraph(inputFile, rotLib, atLib, pairLib);
 
    cout << "get accuracy " << endl;

    for(int i=0;i<keys.size();i++){
        accuracyList.push_back(graph->keyAccuracy(keys[i]));
    }

    ofstream out;
    out.open(outputFile.c_str(), ios::out);
    for(int i=0;i<keys.size();i++){
        out << keys[i] << " " << eneList[i] << " " << accuracyList[i] << endl;
    }
    out.close();

    delete graph;
    delete pairLib;
    delete rotLib;
    delete atLib;

}