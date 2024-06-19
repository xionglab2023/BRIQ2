/**
 * @file bmAssign.cpp
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief Main function of bmAssign, A motif identification tool by community detection of NuGraph
 * @version 0.1
 * @date 2023-12-28
 * 
 * @copyright Copyright (c) 2023
 * 
 */

// #include "briqxmodule/WrapperInfomap.h"
#include "briqxmodule/MotifAssigner.h"
#include "model/AssignRNASS.h"
#include "tools/CmdArgs.h"
#include "tools/ThreadPool.h"
#include <filesystem>
#include <iterator>

using namespace NSPmodel;
using namespace NSPbm;
using namespace NSPtools;
using namespace NSPthread;
using namespace std;

void printHelp() {
    cout << "Usage:" <<endl;
    cout << "bmassign -h" << endl;
    cout << "bmassign -i PDB -o outPath -op outPDBprefix" << endl;
    cout << "bmassign -n NumThreads -f PDBList -o outPath -op outPDBprefix" << endl;
    cout <<
    "If a single PDB file is provided, the assigned motifs are written to outPath/outPDBprefix-x.pdb." << endl;
    cout <<
    "If a list of PDB files are provided, motifs are written to outPath/InputPDBFileName/outPDBprefix-x.pdb." << endl;
}

/**
 * @brief Generate NuGraph Input and run MotifAssigner
 * 
 * @return int 
 */
int main(int argc, char** argv) {
    if(argc == 1) {
        printHelp();
        return EXIT_SUCCESS;
    }
    CmdArgs cmdArgs{argc, argv};
    string listFile, pdbIn;

    filesystem::path outPath{};
    filesystem::path outPDBprefix{"motif"};

    if(cmdArgs.specifiedOption("-h")) {
        printHelp();
        return EXIT_SUCCESS;
    }

    if(cmdArgs.specifiedOption("-o")) {
        outPath = cmdArgs.getValue("-o");
    }
    if(cmdArgs.specifiedOption("-op")) {
        outPDBprefix = cmdArgs.getValue("-op");
    }

    if(cmdArgs.specifiedOption("-i")) {
        if(cmdArgs.specifiedOption("-n")) {
            cout<<"[Warning][Main] OptionIgnored: Option -n is ignored when -i is specified." << endl;
        }
        if(cmdArgs.specifiedOption("-f")) {
            cout<<"[Warning][Main] OptionIgnored: Option -f is ignored when -i is specified." << endl;
        }
        pdbIn = cmdArgs.getValue("-i");
        RNAPDB* rnaPdb = new RNAPDB(pdbIn);
        rnaPdb->DNAToRNA();
        AtomLib* atl = new AtomLib();
        AssignRNASS rss(rnaPdb, atl);
        ofstream tmpInput;
        string tmpInputFileName = outPath.string() + "/tmpIn";
        tmpInput.open(tmpInputFileName, ios::out);
        if (! tmpInput.is_open())
        {
            throw "[Error] Fail to open file " + tmpInputFileName;
        }
        tmpInput << "task analysis" <<endl;
        tmpInput << "pdb " << pdbIn << endl;
        tmpInput << "seq " << rss.seq << endl;
        tmpInput << "sec " << rss.ssSeq << endl;
        string cntStr = "-";
        cntStr.resize(rss.seq.size(),'-');
        int nbreak = rss.breakList.size();
        for(int i=0; i<nbreak; i++) {
            cntStr[rss.breakList[i]] = '|';
        }
        cntStr[rss.seq.size()-1] = '|';
        tmpInput << "cnt " << cntStr << endl;
        string cstStr = ".";
        cstStr.resize(rss.seq.size(),'.');
        tmpInput << "cst " << cstStr << endl;
        tmpInput.close();
        delete rnaPdb;
        RotamerLib* rtl = new RotamerLib();
        BasePairLib* bpl = new BasePairLib();
        RnaEnergyTable* et = new RnaEnergyTable();
        et->loadEnergyWithout6D();
        NuGraph* pNuGragh = new NuGraph(tmpInputFileName, rtl, atl, bpl, et, 1);
        pNuGragh->initInfo->printPDB(
            outPath.string() + "/" + outPDBprefix.string() + "-full" + ".pdb");
        MotifAssigner* mtfa = new MotifAssigner(pNuGragh);
        ofstream outCSV;
        string csvFileName = outPath.string() + "/" + outPDBprefix.string() + "-EdgeWeights" + ".csv";
        outCSV.open(csvFileName, ios::out);
        if(!outCSV.is_open()) {
            throw "[Error] Fail to open file" + csvFileName;
        }
        mtfa->writeEdgeWeight(outCSV);
        outCSV.close();
        
        mtfa->bySeed(outPath.string(), outPDBprefix.string());
        int nMotif = mtfa->motifs.size();
        cout << "find " << nMotif << " motifs" << endl;

        string sepFile = outPath.string() + "/" + outPDBprefix.string() + ".sep";
        ofstream ofSep;
        ofSep.open(sepFile, ios::out);
        if(!ofSep.is_open()) {
            throw "[Error] Can not open file " + sepFile;
        }
        ofSep << "Motif,Sep" << endl;
        for(int i=0; i<nMotif; i++) {
            string motifName = outPDBprefix.stem().string() + "-" + to_string(i);
            mtfa->motifs[i]->setName(motifName);
            mtfa->motifs[i]->writePDB(
                outPath.string() + "/" + outPDBprefix.string() + "-" + to_string(i) + ".pdb");
            ofSep << motifName<<",";
            copy(mtfa->motifs[i]->sepVec.begin(), mtfa->motifs[i]->sepVec.end(), 
                ostream_iterator<int>(ofSep, ","));
            ofSep << endl;
        }
        ofSep.close();

/*  disable redundant seqFrag and IndexInFull output
        string seqFile = outPath.string() + "/" + outPDBprefix.string() + ".seqFrag";
        ofstream ofSeqFrag;
        ofSeqFrag.open(seqFile, ios::out);
        if(!ofSeqFrag.is_open()) {
            throw "[Error] Can not open file " + seqFile;
        }
        mtfa->writeSeqFrag(ofSeqFrag);
        ofSeqFrag.close();

        string ndxFile = outPath.string() + "/" + outPDBprefix.string() + ".ndx";
        ofstream ofNdx;
        ofNdx.open(ndxFile, ios::out);
        if(!ofNdx.is_open()) {
            throw "[Error] Can not open file " + ndxFile;
        }
        mtfa->writeNodeIndexInFull(ofNdx);
        ofNdx.close();
*/

        string bmFile = outPath.string() + "/" + outPDBprefix.string() + ".bm";
        ofstream ofBm;
        ofBm.open(bmFile, ios::out);
        if(!ofBm.is_open()) {
            throw "[Error] Can not open file " + bmFile;
        }
        mtfa->writeMotif(ofBm, rss.ssSeq);
        ofBm.close();

        delete atl;
        delete rtl;
        delete bpl;
        delete pNuGragh;
        delete mtfa;
    }
    return EXIT_SUCCESS;
}
