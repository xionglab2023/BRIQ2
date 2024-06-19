/**
 * @file findRedundantMotif.cpp
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief find redundant motifs from a given list of pdbs
 * @version 0.1
 * @date 2024-06-15
 * 
 * @copyright Copyright (c) 2024 XLAB
 * 
 * @addtogroup BriqxModule
 */

#include "geometry/RMSD.h"
#include "tools/CmdArgs.h"
#include "model/StructureModel.h"
#include <filesystem>

using namespace NSPgeometry;
using namespace NSPmodel;
using namespace NSPtools;
using namespace std;

void printHelp() {
    cout << "Usage:" << endl;
    cout << "findRedundantMotif -h" << endl;
    cout << "findRedundantMotif -f PDBList -o outPath" << endl;
    cout << "findRedundantMotif -f PDBList -o outPath -c rmsdCutoff" << endl;
}

int main(int argc, char** argv) {
    if(argc == 1) {
        printHelp();
        return EXIT_SUCCESS;
    }
    CmdArgs cmdArgs{argc, argv};
    filesystem::path pdbList{"PDBList"}, outPath{"."}; 
    double rmsdC = 1.5;

    if(cmdArgs.specifiedOption("-h")) {
        printHelp();
        return EXIT_SUCCESS;
    }

    if(cmdArgs.specifiedOption("-o")) {
        outPath = cmdArgs.getValue("-o");
    }

    if(cmdArgs.specifiedOption("-f")) {
        pdbList = cmdArgs.getValue("-f");
    }

    if(cmdArgs.specifiedOption("-c")) {
        istringstream ss(cmdArgs.getValue("-c"));
        ss >> rmsdC;
    }

    ifstream input;
    input.open(pdbList, ios::in);
    if (! input.is_open())
    {
        cerr << "[Error][Main] fail to open file " + pdbList.string() + '\n';
        return EXIT_FAILURE;
    }
    string s;
    map<string, filesystem::path>* pdbPathMap = new map<string, filesystem::path>;
    while(getline(input,s)) {
        filesystem::path path0{s};
        pdbPathMap->emplace(path0.stem().string(), path0);
    }
    vector<string*>* pdbids = new vector<string*>;
    pdbids->resize(pdbPathMap->size());
    int i1=0;
    for(auto& it1 : *pdbPathMap) {
        pdbids->at(i1) = new string{it1.first};
        i1++;
    }
    cout << "[Info] Read " << pdbids->size() << " pdbs" << endl;

    vector<string*>* pdbidsNew = new vector<string*>, *pdbidsTmp;
    vector<string*>* pdbidNR = new vector<string*>;  // pdbids non-redundant
    map<string*, vector<string*>*>* redundantMap = new map<string*, vector<string*>*>; 
    while(pdbids->size() > 1) {
        int lp = pdbids->size();
        pdbidsNew->clear();
        string PDBidA = *(pdbids->at(0));
        RNAPDB rnaPDBA(pdbPathMap->at(PDBidA), PDBidA);
        int seqALen = rnaPDBA.getBaseList().size();
        string seqA(seqALen, 'X');
        vector<XYZ> coordsA;
        for(int i=0; i<seqALen; i++) {
            seqA[i] = rnaPDBA.getBaseList()[i]->baseType;
            vector<Atom*>* atomList1 = rnaPDBA.getBaseList()[i]->getAtomList();
            int natom1 = atomList1->size();
            for(int j=0; j<natom1; j++) {
                coordsA.emplace_back(atomList1->at(j)->coord);
            }
        }
        for(int j=1;j<lp;j++) {
            string PDBidB = *(pdbids->at(j));
            RNAPDB rnaPDBB(pdbPathMap->at(PDBidB), PDBidB);
            int seqBLen = rnaPDBB.getBaseList().size();
            if(seqBLen != seqALen) {
                pdbidsNew->emplace_back(pdbids->at(j));
                continue;
            }
            
            string seqB(seqBLen, 'X');
            vector<XYZ> coordsB;
            for(int i=0; i<seqBLen; i++) {
                seqB[i] = rnaPDBB.getBaseList()[i]->baseType;
                vector<Atom*>* atomList1 = rnaPDBB.getBaseList()[i]->getAtomList();
                int natom1 = atomList1->size();
                for(int j=0; j<natom1; j++) {
                    coordsB.emplace_back(atomList1->at(j)->coord);
                }
            }
            if(seqA != seqB) {
                pdbidsNew->emplace_back(pdbids->at(j));
                continue;
            }
            double rmsd1 = rmsd(coordsA, coordsB);
            if(rmsd1 >= rmsdC) {
                pdbidsNew->emplace_back(pdbids->at(j));
                continue;
            }
            try {
                redundantMap->at(pdbids->at(0))->emplace_back(pdbids->at(j));
            } catch(out_of_range) {
                vector<string*>* v = new vector<string*>{pdbids->at(j)};
                redundantMap->emplace(pdbids->at(0), v);
            }    
        }
        pdbidNR->emplace_back(pdbids->at(0));

        pdbidsTmp = pdbids;
        pdbids = pdbidsNew;
        pdbidsNew = pdbidsTmp;
    }
    if(pdbids->size() == 1) pdbidNR->emplace_back(pdbids->at(0));
    delete pdbids;
    delete pdbidsNew;

    //output
    cout << "[Info] find " << pdbidNR->size() << " non-redundant PDBs, in which " <<
        redundantMap->size() << " PDBs have redundant analogs." << endl;
    filesystem::path outNRList, outRMapList;  // output for non-redundant pdbid list, and for redundancy mapping list
    outNRList = outPath / "nonRedundant.list";
    outRMapList = outPath / "redundantMap.csv";

    ofstream output;
    output.open(outNRList.string(), ios::out);
    if (! output.is_open())
    {
        cerr << "[Error][Main] fail to open file " + outNRList.string() + '\n';
        return EXIT_FAILURE;
    }
    int ln = pdbidNR->size();
    for(int i=0; i<ln; i++) {
        output << *(pdbidNR->at(i)) << ' ' << pdbPathMap->at(*(pdbidNR->at(i))).string() << '\n';
    }
    output.close();

    output.open(outRMapList.string(), ios::out);
    if (! output.is_open())
    {
        cerr << "[Error][Main] fail to open file " + outRMapList.string() + '\n';
        return EXIT_FAILURE;
    }
    for(auto& it1 : *redundantMap) {
        output << *(it1.first);
        int ln = it1.second->size();
        for(int i=0; i<ln; i++) {
            output << "," << *(it1.second->at(i));
        }
        output << '\n';
    }
    output.close();
    int l1 = pdbidNR->size();
    for(int i =0; i<l1; i++) {
        delete pdbidNR->at(i);  // strings originally allocated by pdbids now reside in redundantMap and pdbidNR
    }
    delete pdbidNR;
    for(auto& it1 : *redundantMap) {
        int l1 = it1.second->size();
        for(int i =0; i<l1; i++) {
            delete it1.second->at(i);  // strings originally allocated by pdbids now reside in redundantMap and pdbidNR
        }
        delete it1.second;
    }
    delete redundantMap;
}