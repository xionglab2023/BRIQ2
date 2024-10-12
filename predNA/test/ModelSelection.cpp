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
	cout << "briqx_modelSelection -list $PDBLIST -n $MODELNUM -cutoff $RMSDCUTOFF -out $OUTPUT" <<endl;
	cout << "briqx_modelSelection -tre $TREEINFOFILE -n $MODELNUM -cutoff $RMSDCUTOFF -out $OUTPUT" << endl;
}

void selectModel(vector<RNAPDB>& pdbList, int n, double cutoff, const string& output){

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

    string listFile = cmdArgs.getValue("-list");
    int n = atoi(cmdArgs.getValue("-key").c_str());
	double cutoff = atof(cmdArgs.getValue("-cutoff").c_str());
	
	ifstream file, f2;
	file.open(listFile.c_str(), ios::in);
	vector<string> pdbFileList;
	vector<RNAPDB*> pdbList;
	vector<double> eneList;

	string s;
	while(getline(file, s)){
		f2.open(s);
		if(!f2.is_open)
	}


}
