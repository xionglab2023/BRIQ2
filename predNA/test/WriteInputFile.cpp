#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include "model/StructureModel.h"
#include "predNA/NuGraph.h"
#include "tools/CmdArgs.h"
#include "tools/ThreadPool.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredNA;
using namespace std;
using namespace NSPtools;


void printHelp(){
	cout << "Usage: " << endl;
	cout << "briqx_keyInput -in $OLDINPUTFILE -key $KEYFILE -out $OUTPUT" <<endl;
}

int main(int argc, char** argv){

    
    if(argc == 1) {
        printHelp();
        return EXIT_SUCCESS;
    }

    CmdArgs cmdArgs{argc, argv};
	if(cmdArgs.specifiedOption("-h")){
		printHelp();
        return EXIT_SUCCESS;
	}

    string oldInput = cmdArgs.getValue("-in");
    string keyList = cmdArgs.getValue("-key");
    string output = cmdArgs.getValue("-out");

    vector<string> keys;
    ifstream file;
    file.open(keyList.c_str(), ios::in);
    string s;
    vector<string> spt;
    while(getline(file, s)){
        splitString(s, " ", &spt);
        keys.push_back(spt[0]);
    }
    file.close();
    
    NSPtools::InputParser input(oldInput);

    ofstream out;
    char xx[200];
    for(int i=0;i<keys.size();i++){
        sprintf(xx, "%s-%d", output.c_str(), i);
        string outfile = string(xx);
        out.open(outfile.c_str(), ios::out);
        out << "task predict" << endl;

        if(input.specifiedOption("pdb"))
            out << "pdb " << input.getValue("pdb") << endl;
        if(input.specifiedOption("seq"))
            out << "seq " << input.getValue("seq") << endl;
        if(input.specifiedOption("sec"))
            out << "sec " << input.getValue("sec") << endl;
        if(input.specifiedOption("cnt"))
            out << "cnt " << input.getValue("cnt") << endl;
        if(input.specifiedOption("cst"))
        
            out << "cst " << input.getValue("cst") << endl;

        out << "key " << keys[i] << endl;
        out.close();
    }

}