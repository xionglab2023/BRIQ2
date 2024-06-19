/**
 * @file ExtractValidRNAPDB.cpp
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief Extract valid RNA bases from RCSB PDB files
 * @version 0.1
 * @date 2024-04-30
 * 
 * @copyright Copyright (c) 2024 XLAB
 * 
 */

#include "model/StructureModel.h"
#include "tools/CmdArgs.h"
#include <filesystem>

using namespace NSPmodel;
using namespace NSPtools;
using namespace std;

void printHelp() {
    cout <<"Usage:" << endl;
    cout <<"extractValidRNAPDB -h" << endl;
    cout <<"extractValidRNAPDB -i inputPDB -o outputPDB" << endl;
}

int main(int argc, char** argv) {
    if(argc == 1) {
        printHelp();
        return EXIT_SUCCESS;
    }
    CmdArgs cmdArgs{argc, argv};
    if(cmdArgs.specifiedOption("-h")) {
        printHelp();
        return EXIT_SUCCESS;
    }

    filesystem::path outPDB{};
    filesystem::path inPDB{};

    if(cmdArgs.specifiedOption("-i")) {
        inPDB = cmdArgs.getValue("-i");
    } else {
        cerr <<"inPutPDB not specified" << endl;
        return EXIT_FAILURE;
    }
    if(cmdArgs.specifiedOption("-o")) {
        outPDB = cmdArgs.getValue("-o");
    } else {
        cerr <<"outPDB not specified" << endl;
        return EXIT_FAILURE;
    }

    RNAPDB* rnaPDB = new RNAPDB(inPDB.string());
    AtomLib* atl = new AtomLib();

    auto vbl = rnaPDB->getValidBaseList(atl);
    
    RNAPDB* validPDB = new RNAPDB();
    for(RNABase* base1 : vbl) {
        if(validPDB->getChain(base1->chainID)==NULL) {
            RNAChain* curChain = new RNAChain(base1->chainID);
            validPDB->getChains().emplace_back(curChain);
        }
        validPDB->getChain(base1->chainID)->addBase(base1);
    }

    ofstream outs;
    outs.open(outPDB, ios::out);
    if (! outs.is_open())
    {
        cerr << "[Error] Fail to open file " << outPDB.string() << endl;
        return EXIT_FAILURE;
    }
    validPDB->printPDBFormat(outs);
    cout << "Valid seqLen: " << vbl.size() << endl;
    outs.close();

    delete rnaPDB;
    delete atl;
    delete validPDB;
    return EXIT_SUCCESS;
}