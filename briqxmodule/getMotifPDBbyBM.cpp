/**
 * @file getMotifPDBbyBM.cpp
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief  Extract motif structures from a full structure PDB accirding to given bm
 * @version 0.1
 * @date 2024-06-06
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "model/StructureModel.h"
#include "tools/CmdArgs.h"
#include "tools/InputParser.h"
#include <filesystem>
#include <algorithm>
// #define DEBUG

using namespace NSPmodel;
using namespace NSPtools;
using namespace std;

void printHelp() {
    cout << "Usage:" << endl;
    cout << "getMotifPDB -h" << endl;
    cout << "getMotifPDB -s fullStructurePDB -b bm -o outPDBPath -m motifBasename" << endl;
    cout << "Extract motif PDBs from **fullStrucurePDB** according to the provided **bm**." << endl;
    cout << "All motifs present in **bm** will be extracted. Output PDBs are named after **motifBasename**, " << "\n" << 
    "numbered after motif mask labels in **bm** and written to directory specified by **outPDBPath**." << endl;
}

int main(int argc, char** argv) {
    if(argc == 1) {
        printHelp();
        return EXIT_SUCCESS;
    }
    CmdArgs cmdArgs{argc, argv};

    filesystem::path bmfile{"bm"};
    filesystem::path fullPDB{"fullPDB"};
    filesystem::path outPath{"."};
    string motifBasename = "motif";

    if(cmdArgs.specifiedOption("-h")) {
        printHelp();
        return EXIT_SUCCESS;
    }

    if(cmdArgs.specifiedOption("-o")) {
        outPath = cmdArgs.getValue("-o");
    }
    if(cmdArgs.specifiedOption("-s")) {
        fullPDB = cmdArgs.getValue("-s");
    }
    if(cmdArgs.specifiedOption("-b")) {
        bmfile = cmdArgs.getValue("-b");
    }
    if(cmdArgs.specifiedOption("-m")) {
        motifBasename = cmdArgs.getValue("-m");
    }

    // read bmfile
    InputParser input(bmfile.string());
    input.printOptions();

    string seq = input.getValue("seq");
    string dbn = input.getValue("sec");
    string cnt = input.getValue("cnt");
    if(seq.size() != dbn.size()) {
        cerr << "Inconsistent length between seq and sec" << endl;
        return EXIT_FAILURE;
    }
    if(seq.size() != cnt.size()) {
        cerr << "Inconsistent length between seq and cnt" << endl;
        return EXIT_FAILURE;
    }

    vector<string> tags = input.listOptionTags();
    sort(tags.begin(), tags.end());
    const auto dup = adjacent_find(tags.begin(), tags.end());
    if (dup != tags.end()) {
        cerr << "[Error] duplicate tag: " << *dup << " in bm file" << endl;
        return EXIT_FAILURE;
    }
    map<string,string> motif_mask;
    int ntag = tags.size();
    for(int i=0; i<ntag; i++) {
        if(tags[i]=="seq") continue;
        if(tags[i]=="sec") continue;
        if(tags[i]=="cnt") continue;
        string key = motifBasename + '-' + tags[i].substr(1);
        string mask = input.getValue(tags[i]);
        if(seq.size() != mask.size()) {
            cerr << "Inconsistent length between seq and mask: " << key << endl;
            return EXIT_FAILURE;
        }
        motif_mask.emplace(key, mask);
    }
    // end read bmfile

    // read fullPDB
    RNAPDB* rnaPdb = new RNAPDB(fullPDB.string());
    auto& fullBaseList = rnaPdb->getBaseList();
    int nBase = fullBaseList.size();
    #ifdef DEBUG
        cout << seq << endl;
        cout << "seq size is " << seq.size() << endl;
        cout << "seq length is " << seq.length() << endl;
        // string teststr = "12345";
        // cout << "12345 length is" << teststr.length() << endl;
        // cout << "seq[25]: " << int(seq[25]) << endl;
        // cout << "seq[26]: " << int(seq[26]) << endl;
    #endif
    if (seq.size() > nBase) {
        cerr << "[Error] " << motifBasename <<" full PDB shorter than seq" << endl;
        return EXIT_FAILURE;
    }
    for (auto& it1 : motif_mask) {
        RNAPDB* motifPDB = new RNAPDB();
        for (int i=0; i< seq.size(); i++) {
            if(it1.second[i] == '-') continue;
            RNABase* base1 = fullBaseList[i];
            if(base1->baseType != seq[i]) {
                cerr << "[Error] " << motifBasename<<" fullPDB base " << base1->baseID << base1->baseType << " does not match seq " << seq[i] << endl;
                return EXIT_FAILURE;
            }
            auto curChain = motifPDB->getChain(base1->chainID);
            if(curChain==NULL) {
                RNAChain *newChain = new RNAChain(base1->chainID);
                motifPDB->getChains().emplace_back(newChain);
                newChain->addBase(base1);
            } else {
                curChain->addBase(base1);
            }                
        }

        ofstream outs;
        auto outPDB = outPath / (it1.first + ".pdb");
        outs.open(outPDB.string(), ios::out);
        if (! outs.is_open())
        {
            cerr << "[Error] Fail to open file " << outPDB.string() << endl;
            return EXIT_FAILURE;
        }
        motifPDB->printPDBFormat(outs);
        outs.close();

        delete motifPDB;
    }
    delete rnaPdb;
    return EXIT_SUCCESS;
}