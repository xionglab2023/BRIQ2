/**
 * @file Main.cpp
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief  Main function implementation
 * @version udef
 * @date 2023/08/21
 * 
 * @copyright Copyright (c) 2023 XLAB
 * 
 * modification history :
 * Date:      Version:    Author:
 * Changes:
 */

/**
 * @addtogroup BriqxModule
 * @brief  Main process of BriqxModule alignment
 * 
 * Main process for BriqxModule alignment, inlcuding user interaction and output formating.
 * 
 * @{
 */

#include "briqxmodule/BriqxModule.h"
#include "briqxmodule/BMAlign.h"
#include "briqxmodule/BMScore.h"
#include "tools/CmdArgs.h"
#include "tools/ThreadPool.h"
#include <filesystem>

#define DDMCUTOFF 0.7
#define ITERSTOPAT -1e-6

using namespace NSPtools;
using namespace NSPbm;
using namespace NSPthread;
using namespace std;

/**
 * @brief Workhorse function to align and calculate alignment score for two BriqxModules, single thread 
 * 
 * @param [in] BMaPDB: PDB file for BriqxModule a, required.
 * @param [in] BMbPDB: PDB file for BriqxModule b, required.
 * @param [in] bpl: BasePairLib object for basepair type and energy determination, required.
 * @param [in] atl: AtomLib object for BasePair determination, required.
 * @param [out] score: Score of the best alignment, required. 
 * @param [out] alignVec: Vector of the aligned bases in the best alignment, required. 
 * @param [out] NormedSc: Normalized Score of the best alignment, required. 
 * @param [in] BMaSel: Select string from the PDB for BMa, optional.
 * @param [in] BMbSel: Select string from the PDB for BMb, optional.
 * @param [in] BMaName: Name of BMa, optional.
 * @param [in] BMbName: Name of BMb, optional.
 * @param [in] outPath: Path to write output files, optional.
 * @param [in] outfileAln: Name of output pairwise alignment file, no suffix, in a2m(fastsa) format, optional.
 * @param [in] outfilePDB: Name of output aligned BMb coordinates, no suffix, in PDB format, optional.
 * 
 * @return success code.
 */
int alignAndScore(const string& BMaPDB, const string& BMbPDB, BasePairLib& bpl, AtomLib& atl,
    double& score, vector<array<RNABase*, 2> >& alignVec, double& normedSc,
    const BMSelection& BMaSel = 0, const BMSelection& BMbSel = 0,
    const string& BMaName = "BMa", const string& BMbName = "BMb",
    const string& outPath = NULL,
    const string& outfileAln = "align.aln", const string& outfilePDB = "alignedB.pdb")
{
    
    const double DDMCutoff = DDMCUTOFF;
    const double iterStopAt = ITERSTOPAT;

    BriqxModule BMa(BMaPDB, bpl, atl, BMaName, BMaSel);
    BriqxModule BMb(BMbPDB, bpl, atl, BMbName, BMbSel);
    #ifdef DEBUG
    BMa.printPairEneMap(cout);
    BMb.printPairEneMap(cout);  
    #endif

    if(BMa.getNStrand() != BMb.getNStrand()) {
        string outstr = "[Result] nStrand Unmatch: " +\
            BMaName + ": " + to_string(BMa.getNStrand()) + ", " + \
            BMbName + ": " + to_string(BMb.getNStrand()) + "\n";
        cout << outstr;
        score = -1;
        normedSc = -1;
        alignVec.clear();
        return EXIT_SUCCESS;
    }

    map<array<BasePair*, 2>, double> DDMMatrix{};

    BMa.calcDDMMatrix(BMb, DDMMatrix);

    BMScore scer(BMa, BMb, DDMMatrix);
    #ifdef DEBUG
        scer.printBasePairWeight(cout);
    #endif
    double idealSc = scer.idealScore();

    vector<pair<array<BasePair*, 2>, double> > sortedDDM{};

    // scer.sortBppByScore(sortedDDM);  // 应该以DDM排序，因为需要遍历所有DDM<DDMcutoff的 BasePair 对
    utils::sortMapByValue(DDMMatrix, sortedDDM, false);

    int i = 0;
    auto* p1 = &sortedDDM[i];
    double bestScore = 0;
    BMAlign<vector<array<RNABase*,2> > >* bestAlign = nullptr;
    BMAlign<vector<array<BasePair*,2> > >* bestAlignBP = nullptr;
    // while(DDMMatrix.at(p1->first) <= DDMCutoff) {
    while(p1->second <= DDMCutoff) {
        #ifdef DEBUG
            string outstr = "[DEBUG][Info] Evaluating initial alignment by pair " + to_string(i) + ": " + \
                (p1->first)[0]->print() + " To " + (p1->first)[1]->print() + "\n";
            cout << outstr;
            int itCount = 0;
        #endif //DEBUG
        auto* iniAlign = new BMAlign<vector<array<BasePair*,2> > >(BMa, BMb, vector<array<BasePair*,2> >{p1->first});
        // #ifdef DEBUG
            // ofstream debugOut;
            // debugOut.open("iniTransformed.pdb", ios::out);
            // auto BMnew = BriqxModule(BMb, iniAlign->getRotTrans(), iniAlign->getAcog(), iniAlign->getBcog(),
            //      bpl, atl, "", true);
            // BMnew.printPDBFormat(debugOut);
            // debugOut.close();
            // BMnew.deepClear();
        // #endif //DEBUG
        if(iniAlign->getAlignVec().size() == 0) {
            #ifdef DEBUG
                cout << "[DEBUG][Info] Stop coz initial alignment get empty alignVec" <<endl;
            #endif //DEBUG
            break;
        }
        BMAlign<vector<array<RNABase*,2> > >* curAlign = nullptr, *lastAlign = nullptr, *nextAlign;
        double curSc = 0, lastSc;
        double nextSc = scer.score(iniAlign->getAlignVec());
        if(curSc-nextSc < iterStopAt) {
            #ifdef DEBUG
                itCount++;
                cout << "[DEBUG][Info] IniBPP: " << i << ", Iteration " << itCount << endl;
            #endif // DEBUG
            lastSc = curSc;
            curSc = nextSc;
            nextAlign = new BMAlign<vector<array<RNABase*,2> > >(BMa, BMb, iniAlign->getAlignVec());
            nextSc = scer.score(nextAlign->getAlignVec());
            #ifdef DEBUG
                outstr = "[DEBUG][Info] New Score: " + to_string(nextSc) + "\n";
                cout << outstr;
            #endif
        }
        while(curSc-nextSc < iterStopAt) {
            #ifdef DEBUG
                itCount++;
                cout << "[DEBUG][Info] IniBPP: " << i << ", Iteration " << itCount << endl;
            #endif // DEBUG
            lastSc = curSc;
            curSc = nextSc;
            lastAlign = curAlign;
            curAlign = nextAlign;
            nextAlign = new BMAlign<vector<array<RNABase*,2> > >(BMa, BMb, curAlign->getAlignVec());
            nextSc = scer.score(nextAlign->getAlignVec());
            #ifdef DEBUG
                outstr = "[DEBUG][Info] New Score: " + to_string(nextSc) + "\n";
                cout << outstr;
            #endif
            if(lastAlign) delete lastAlign;
        }
        delete nextAlign;
        #ifdef DEBUG
            cout << "[DEBUG][Info] Final alignment starting from BPP " << i << ": " << (p1->first)[0]->print() << " To " <<
                    (p1->first)[1]->print() << endl;
            double nSc = curSc/idealSc;
            string outAlnPath = outPath + "/" + BMb.getBMname() + "_to_" + BMa.getBMname() +\
                "_by_" + (p1->first)[0]->print() + '-' + (p1->first)[1]->print() + ".aln";
            ofstream dbg_output;
            dbg_output.open(outAlnPath, ios::out);
            if (! dbg_output.is_open())
            {
                throw "[Error] Fail to open file " + outAlnPath;
            }
            if(curAlign) {
                curAlign->writeAlignment(dbg_output, curSc, nSc, true);
            } else {
                iniAlign->writeAlignment(dbg_output, curSc, nSc, true);
            }
            dbg_output.close();

            // Writing intermediate aligned conformation
            string outPDBPath = outPath + "/" + BMb.getBMname() + "_to_" + BMa.getBMname() +\
                "_by_" + (p1->first)[0]->print() + '-' + (p1->first)[1]->print() + ".pdb";
            cout<<"[DEBUG][Info] Write transformed "<<BMb.getBMname()<<" to " << outPDBPath <<endl;
            dbg_output.open(outPDBPath, ios::out);
            if (! dbg_output.is_open())
            {
                throw "[Error] Fail to open file " + outPDBPath;
            }
            BriqxModule* BMdbg;
            if(curAlign) {
                BMdbg = new BriqxModule(BMb, curAlign->getRotTrans(), curAlign->getAcog(), curAlign->getBcog(),
                bpl, atl, "", true);
            } else {
                BMdbg = new BriqxModule(BMb, iniAlign->getRotTrans(), iniAlign->getAcog(), iniAlign->getBcog(),
                bpl, atl, "", true);
            }
            BMdbg->printPDBFormat(dbg_output);
            dbg_output.close();
            BMdbg->deepClear();
            delete BMdbg;
        #endif
        if(curSc > bestScore) {
            bestScore = curSc;
            #ifdef DEBUG
                outstr = "[DEBUG][Info] Refreshing, best score now from iniBPP " + to_string(i) + " with value " +\
                    to_string(bestScore) + "\n";
                cout << outstr;
            #endif //DEBUG
            if(curAlign) {
                if(bestAlign) delete bestAlign;
                bestAlign = curAlign;
            } else {
                if(bestAlignBP) delete bestAlignBP;
                bestAlignBP = iniAlign;
            }
        } else {
            delete iniAlign;
            delete curAlign;
        }

        try {   
            p1 = &sortedDDM.at(++i);
        } catch(out_of_range) {
            break;
        }
    }
    
    if(bestScore == 0) {
        cout << "[Result] NoAlignment: " << BMaName << "-" << BMbName << endl;
        score = -2;
        normedSc = -2;
        alignVec.clear();
        return EXIT_SUCCESS;
    }
    score = bestScore;
    alignVec = bestAlign? bestAlign->getAlignVec(): bestAlignBP->getAlignVec();
    normedSc = score/idealSc;
 
    cout<<"[Result] "<< BMaName <<"-"<<BMbName<<" Base-Base alignment completed" <<endl;
    cout<<"[Result] Alignment score:" << to_string(score) << endl;
    cout<<"[Result] Normed score:" << to_string(normedSc) << endl;

    if(outPath != "") {
        filesystem::path pthAln(outfileAln);
        auto pAStem = pthAln.stem();
        auto pAExt = pthAln.extension();
        cout<<"[Info] Write alignment on "<<BMaName<<" to " <<outPath<<"/"<<pAStem.string()<<"_onA"<<pAExt.string()<<
            endl;
        ofstream output;
        string filePath = outPath+"/"+outfileAln+"_onA";
        output.open(filePath, ios::out);
        if (! output.is_open())
        {
            throw "[Error] Fail to open file " + filePath;
        }
        if(bestAlign) {
            bestAlign->writeAlignment(output, score, normedSc, true);
        } else {
            bestAlignBP->writeAlignment(output, score, normedSc, true);
        }
        output.close();
        cout<<"[Info] Write alignment on "<<BMaName<<" to " <<outPath<<"/"<<pAStem.string()<<"_onB."<<
            pAExt.string()<<endl;
        filePath = outPath+"/"+outfileAln+"_onB";
        output.open(filePath, ios::out);
        if (! output.is_open())
        {
            throw "[Error] Fail to open file " + filePath;
        }
        if(bestAlign) {
            bestAlign->writeAlignment(output, score, normedSc, false);
        } else {
            bestAlignBP->writeAlignment(output, score, normedSc, false);
        }
        output.close();
        cout<<"[Info] Write transformed "<<BMb.getBMname()<<" to " <<outPath<<"/"<<outfilePDB<<endl;
        filePath = outPath + "/" + outfilePDB;
        output.open(filePath, ios::out);
        if (! output.is_open())
        {
            throw "[Error] Fail to open file " + filePath;
        }
        BriqxModule* BMnew;
        if(bestAlign) {
            BMnew =  new BriqxModule(BMb, bestAlign->getRotTrans(), bestAlign->getAcog(), bestAlign->getBcog(),
            bpl, atl, "", true);
        } else {
            BMnew = new BriqxModule(BMb, bestAlignBP->getRotTrans(), bestAlignBP->getAcog(), bestAlignBP->getBcog(),
            bpl, atl, "", true);
        }
        BMnew->printPDBFormat(output);
        output.close();
        BMnew->deepClear();
        delete BMnew;
    }

    BMa.deepClear();
    BMb.deepClear();
    if(bestAlign) {
        delete bestAlign;
    } else {
        delete bestAlignBP;
    }
    return EXIT_SUCCESS;
}


/**
 * @brief Workhorse function to align and calculate alignment score for two BriqxModules, multi-thread 
 * 
 * @param [in] BMaPDB: PDB file for BriqxModule a, required.
 * @param [in] BMbPDB: PDB file for BriqxModule b, required.
 * @param [in] bpl: BasePairLib object for basepair type and energy determination, required.
 * @param [in] atl: AtomLib object for BasePair determination, required.
 * @param [in] jid: id of current job.
 * @param [out] scoreMap: Score of the best alignment, required. In multi-thread version, it is a pre-allocated
 *      map<array<string, 2>, double> where the strings are BMaName and BMbName
 * @param [out] normedScMap: Normalized Score of the best alignment, required. 
 * @param [in] BMaSel: Select string from the PDB for BMa, optional.
 * @param [in] BMbSel: Select string from the PDB for BMb, optional.
 * @param [in] BMaName: Name of BMa, optional.
 * @param [in] BMbName: Name of BMb, optional.
 * @param [in] outPath: Path to write output files, optional.
 * @param [in] outfileAln: Name of output pairwise alignment file, no suffix, in a2m(fastsa) format, optional.
 * @param [in] outfilePDB: Name of output aligned BMb coordinates, no suffix, in PDB format, optional.
 * 
 * @return success code.
 */
int alignAndScoreMT(string BMaPDB, string BMbPDB, shared_ptr<BasePairLib> bpl, shared_ptr<AtomLib> atl, size_t jid,
    shared_ptr<map<array<string, 2>, double> > scoreMap, shared_ptr<map<array<string, 2>, double> > normedScMap,
    shared_ptr<BMSelection> BMaSel = nullptr, shared_ptr<BMSelection> BMbSel = nullptr,
    string BMaName = "BMa", string BMbName = "BMb",
    string outPath = NULL,
    string outfileAln = "align.aln", string outfilePDB = "alignedB.pdb")
{
    
    const double DDMCutoff = DDMCUTOFF;
    const double iterStopAt = ITERSTOPAT;

    double score, normedSc;
    vector<array<RNABase*, 2> > alignVec;

    shared_ptr<BMSelection> emptyBMS(new BMSelection);
    if(BMaSel == nullptr) {
        BMaSel = emptyBMS;
    }
    if(BMbSel == nullptr) {
        BMbSel = emptyBMS;
    }
    BriqxModule BMa(BMaPDB, *bpl, *atl, BMaName, *BMaSel);
    BriqxModule BMb(BMbPDB, *bpl, *atl, BMbName, *BMbSel);

    if(BMa.getNStrand() != BMb.getNStrand()) {
        string outstr = "[Result][" + to_string(jid) + "] nStrand Unmatch: " +\
            BMaName + ": " + to_string(BMa.getNStrand()) + ", " + \
            BMbName + ": " + to_string(BMb.getNStrand()) + "\n";
        cout << outstr;
        array<string,2> alnKey{BMaName, BMbName};
        scoreMap->at(alnKey) = -1;
        normedScMap->at(alnKey) = -1;
        return EXIT_SUCCESS;
    }

    map<array<BasePair*, 2>, double> DDMMatrix{};

    BMa.calcDDMMatrix(BMb, DDMMatrix);

    BMScore scer(BMa, BMb, DDMMatrix);

    vector<pair<array<BasePair*, 2>, double> > sortedDDM{};

    // scer.sortBppByScore(sortedDDM);
    utils::sortMapByValue(DDMMatrix, sortedDDM, false);

    int i = 0;
    auto* p1 = &sortedDDM[i];
    double bestScore = 0;
    BMAlign<vector<array<RNABase*,2> > >* bestAlign = nullptr;
    BMAlign<vector<array<BasePair*,2> > >* bestAlignBP = nullptr;
    // while(DDMMatrix.at(p1->first) <= DDMCutoff) {
    while(p1->second <= DDMCutoff) {
        auto* iniAlign = new BMAlign<vector<array<BasePair*,2> > >(BMa, BMb, vector<array<BasePair*,2> >{p1->first});
        #ifdef DEBUG
            string outstr = "[DEBUG][Info][" + to_string(jid) + "] Evaluating initial alignment by pair " +\
                            to_string(i) + "\n";
            cout << outstr ;
            int itCount = 0;
            // ofstream debugOut;
            // debugOut.open("iniTransformed.pdb", ios::out);
            // auto BMnew = BriqxModule(BMb, iniAlign->getRotTrans(), iniAlign->getAcog(), iniAlign->getBcog(),
            //      bpl, atl, "", true);
            // BMnew.printPDBFormat(debugOut);
            // debugOut.close();
            // BMnew.deepClear();
        #endif //DEBUG
        if(iniAlign->getAlignVec().size() == 0) {
            #ifdef DEBUG
                outstr = "[DEBUG][Info][" + to_string(jid) + "] Stop coz initial alignment get empty alignVec\n";
                cout << outstr; // << operator is not thread-safe, i.e. we can use it atomically only
            #endif //DEBUG
            delete iniAlign;
            break;
        }
        BMAlign<vector<array<RNABase*,2> > >* curAlign = nullptr, *lastAlign = nullptr, *nextAlign;
        double curSc = 0, lastSc;
        double nextSc = scer.score(iniAlign->getAlignVec());
        if(curSc-nextSc < iterStopAt) {
            #ifdef DEBUG
                itCount++;
                outstr = "[DEBUG][Info][" + to_string(jid) + "] IniBP: " + to_string(i) + ", Iteration " + \
                         to_string(itCount) + "\n";
                cout << outstr;
            #endif // DEBUG
            lastSc = curSc;
            curSc = nextSc;
            nextAlign = new BMAlign<vector<array<RNABase*,2> > >(BMa, BMb, iniAlign->getAlignVec());
            nextSc = scer.score(nextAlign->getAlignVec());
        }
        while(curSc-nextSc < iterStopAt) {
            #ifdef DEBUG
                itCount++;
                outstr = "[DEBUG][Info][" + to_string(jid) + "] IniBP: " + to_string(i) + ", Iteration " + \
                to_string(itCount) + "\n";
                cout << outstr;
            #endif // DEBUG
            lastSc = curSc;
            curSc = nextSc;
            lastAlign = curAlign;
            curAlign = nextAlign;
            nextAlign = new BMAlign<vector<array<RNABase*,2> > >(BMa, BMb, curAlign->getAlignVec());
            nextSc = scer.score(nextAlign->getAlignVec());
            if(lastAlign) delete lastAlign;
        }
        delete nextAlign;
        if(curSc > bestScore) {
            bestScore = curSc;
            #ifdef DEBUG
                outstr = "[DEBUG][Info][" + to_string(jid) + "] Refreshing, best score now from iniBP " + \
                         to_string(i) + "\n";
                cout << outstr; 
            #endif //DEBUG
            if(curAlign) {
                if(bestAlign) delete bestAlign;
                delete iniAlign;
                bestAlign = curAlign;
            } else {
                if(bestAlignBP) delete bestAlignBP;
                bestAlignBP = iniAlign;
            }
        } else {
            delete iniAlign;
            delete curAlign;
        }

        try {   
            p1 = &sortedDDM.at(++i);
        } catch(out_of_range) {
            break;
        }
    }
    
    if(bestScore == 0) {
        string outstr = "[Result][" + to_string(jid) + "] NoAlignmnet: " + BMaName + " vs " + BMbName + "\n";
        cout << outstr;
        score = -2;
        normedSc = -2;
        array<string,2> alnKey{BMaName, BMbName};
        scoreMap->at(alnKey) = score;
        normedScMap->at(alnKey) = normedSc;
        alignVec.clear();
        return EXIT_SUCCESS;
    }
    score = bestScore;
    alignVec = bestAlign? bestAlign->getAlignVec(): bestAlignBP->getAlignVec();
    double idealSc = scer.idealScore();
    normedSc = score/idealSc;

    string outstr = "[Result][" + to_string(jid) + "] " + BMa.getBMname() + "-" + BMb.getBMname() +\
        " Base-Base alignment completed\n";
    cout<<outstr;
    outstr = "[Result][" + to_string(jid) + "] Alignment score:" + to_string(score) + "\n";
    cout<<outstr;
    outstr = "[Result][" + to_string(jid) + "] Normed score:" + to_string(normedSc) + "\n";
    cout<<outstr;
    array<string,2> alnKey{BMaName, BMbName};
    scoreMap->at(alnKey) = score;
    normedScMap->at(alnKey) = normedSc;

    if(outPath != "") {
        if(outfileAln.compare("//") != 0) {
            filesystem::path pthAln(outfileAln);
            auto pAStem = pthAln.stem();
            auto pAExt = pthAln.extension();
            string outstr = "[Info][" + to_string(jid) + "] Write alignment on " + BMa.getBMname() + " to "  + outPath +\
            "/" + pAStem.string() + "_" + BMaName + "_with_" + BMbName + pAExt.string() + "\n";
            cout << outstr;
            ofstream output;
            string filePath = outPath+"/" + pAStem.string() + "_" + BMaName + "_with_" + BMbName + pAExt.string();
            output.open(filePath, ios::out);
            if (! output.is_open())
            {
                throw "[Error][" + to_string(jid) + "] Fail to open file " + filePath;
            }
            if(bestAlign) {
                bestAlign->writeAlignment(output, score, normedSc, true);
            } else {
                bestAlignBP->writeAlignment(output, score, normedSc, true);
            }
            output.close();
        }
        if(outfilePDB.compare("//") != 0) {
            filesystem::path pthPDB(outfilePDB);
            auto pPStem = pthPDB.stem();
            auto pPExt = pthPDB.extension();
            outstr = "[Info][" + to_string(jid) + "] Write transformed " + BMbName + " to "  + BMaName + " to " + outPath \
                + "/" + pPStem.string() + "_" + BMaName + "_To_" + BMbName + pPExt.string() + "\n";
            cout<<outstr;
            string filePath = outPath + "/" + pPStem.string() + "_" + BMbName + "_To_" + BMaName + pPExt.string();
            ofstream output;
            output.open(filePath, ios::out);
            if (! output.is_open())
            {
                throw "[Error][" + to_string(jid) + "] Fail to open file " + filePath;
            }
            BriqxModule* BMnew;
            if(bestAlign) {
                BMnew =  new BriqxModule(BMb, bestAlign->getRotTrans(), bestAlign->getAcog(), bestAlign->getBcog(),
                *bpl, *atl, "", true);
            } else {
                BMnew = new BriqxModule(BMb, bestAlignBP->getRotTrans(), bestAlignBP->getAcog(), bestAlignBP->getBcog(),
                *bpl, *atl, "", true);
            }
            BMnew->printPDBFormat(output);
            output.close();
            BMnew->deepClear();
            delete BMnew;
        }
    }

    BMa.deepClear();
    BMb.deepClear();
    if(bestAlign) {
        delete bestAlign;
    } else {
        delete bestAlignBP;
    }
    return EXIT_SUCCESS;
}

void printHelp() {
    cout << "Usage:" <<endl;
    cout << "bmalign -a PDBA -b PDBB -o outPath -op outPDBFileName -oa outAlignmentFileName" << endl;
    cout << "bmalign -n NumThreads -f PDBList -o outPath -op outPDBFileName -oa outAlignmentFileName"<<
            " -os outScoreFileName" << endl;
    cout << "bmalign -a PDBA -f PDBList -n NumThreads -o outPath -op outPDBFileName -oa outAlignmentFileName"<<
            " -os outScoreFileName" << endl;
}

/**
 * @brief  Main function of the Program.
 * 
 * @param Options accepted:
 * @param -h: print help
 * @param -a: path to module A PDB file
 * @param -b: Path to module B PDB file 
 * @param -o: Output path
 * @param -op: filename (no path included) of output PDB file recording aligned module B.
 *             default: alignedB.pdb without -f, "//"(suppressed) when -f is specified.
 * @param -oa: filename (no path included) of output alignment files (based on module A and module B, respectively).
 *             default: align.aln without -f, "//"(suppressed) when -f is specified.
 * @param -f: Read input PDB list from given file.
 * @param -os: filename (path included) to write overall scores in csv format, suffixes labeling nStrand and ".csv" will
 * be appended. Default is path from -o with filename "overallScores".
 * @param -n: number of threads to use.
 * @param Options not implemented:
 * @param -d: Read input PDB from given directory.
 * @param -selA: Apply base-based selections from input PDB of module A. 
 * @param -selB: Apply base-based selections from input PDB of module B. 
 * 
 * @return int 
 */
int main(int argc, char** argv) {
    if(argc == 1) {
        printHelp();
        return EXIT_SUCCESS;
    }
    CmdArgs cmdArgs{argc, argv};
    string listFile, pdbA, pdbB, sel;
    filesystem::path outPath{};
    string outfileAln = "align.aln";
    string outfilePDB = "alignedB.pdb";

    double score, normedSc;
    vector<array<RNABase*,2> > alignVec;

    if(cmdArgs.specifiedOption("-h")) {
        printHelp();
        return EXIT_SUCCESS;
    }

    if(cmdArgs.specifiedOption("-o")) {
        outPath = cmdArgs.getValue("-o");
    }
    if(cmdArgs.specifiedOption("-oa")) {
        outfileAln = cmdArgs.getValue("-oa");
    }
    if(cmdArgs.specifiedOption("-op")) {
        outfilePDB = cmdArgs.getValue("-op");
    }
    filesystem::path outfileSc = outPath / "overallScores";
    if(cmdArgs.specifiedOption("-os")) {
        outfileSc = cmdArgs.getValue("-os");
    }
    if(cmdArgs.specifiedOption("-sel")) {
        sel = cmdArgs.getValue("-sel");   // parser for sel string has not been implemented.
    }
    if(cmdArgs.specifiedOption("-a")) {
        pdbA = cmdArgs.getValue("-a");  // -a and -b must be specified at the same time
        if(cmdArgs.specifiedOption("-d")) {
            cout<<"[Warning][Main] OptionIgnored: Option -d is ignored at the presence of -a,-b,-f"<<endl;
        }
        if(cmdArgs.specifiedOption("-b")) {
            pdbB = cmdArgs.getValue("-b");
            BasePairLib bpl;
            AtomLib atl;
            alignAndScore(pdbA, pdbB, bpl, atl, score, alignVec, normedSc, 0, 0,  // for now use NULL(i.e. 0) BMSelection
                "BMa", "BMb", outPath, outfileAln, outfilePDB);
            return EXIT_SUCCESS;
        } else if(cmdArgs.specifiedOption("-f")) {
            listFile = cmdArgs.getValue("-f");  // read PDB file names from listfile
            cout << "[Info][Main] Align Motifs from " + listFile + " to " + pdbA << endl;
            filesystem::path pathA{pdbA};
            string stemA = pathA.stem();
            if(! cmdArgs.specifiedOption("-oa")) {
                outfileAln = "//";
            }
            if(! cmdArgs.specifiedOption("-op")) {
                outfilePDB = "//";
            }
            int nt = 1;
            if(cmdArgs.specifiedOption("-n")) {
                istringstream ss(cmdArgs.getValue("-n"));
                ss >> nt;
            }

            //batch processing
            ifstream input;
            input.open(listFile, ios::in);
            if (! input.is_open())
            {
                throw "[Error][Main] fail to open file " + listFile;
            }
            string s;
            map<int,vector<string> > strandNum2MotifMap;
            shared_ptr<BasePairLib> pbpl(new BasePairLib);
            shared_ptr<AtomLib> patl(new AtomLib);
            while(getline(input,s)) {
                BriqxModule BM0(s, *pbpl, *patl, "BM0", 0, true);
                int nStrand = BM0.getNStrand();
                if(strandNum2MotifMap.contains(nStrand)) {
                    strandNum2MotifMap[nStrand].emplace_back(s);
                } else{
                    vector<string> vec0{s};
                    strandNum2MotifMap.emplace(nStrand, vec0);
                }
            }
            BriqxModule BM0(pdbA, *pbpl, *patl, "BM0", 0, true);
            int nStrandA = BM0.getNStrand();
            input.close();
            cout<<"[Info][Main] List of Motifs Summary:"<<endl;
            cout<<"     nStrand  Population"<<endl;
            for(auto &iter: strandNum2MotifMap) {
                cout<<setw(12)<<iter.first<<setw(12)<<iter.second.size()<<endl;
            }
            shared_ptr<ThreadPool> thrPool(new ThreadPool(nt));
            size_t jid = 0;
            auto strandVec = strandNum2MotifMap.at(nStrandA);
            shared_ptr<map<array<string, 2>, double> > scoreMap(new map<array<string, 2>, double>);
            shared_ptr<map<array<string, 2>, double> > normedScMap(new map<array<string, 2>, double>);
            shared_ptr<BMSelection> emptyBMS(new BMSelection);
            vector<string> stemVec;  // vector of filename w/o suffix
            int lv = strandVec.size();
            string BMaName = stemA;
            for(int i=0;i<lv;i++) {
                filesystem::path path0{strandVec[i]};
                stemVec.emplace_back(path0.stem());
                string BMbName = stemVec[i];
                scoreMap->emplace(array<string,2>{BMaName,BMbName}, 0);
                normedScMap->emplace(array<string,2>{BMaName,BMbName}, 0);
            }
            for(int j=0;j<lv;j++) {
                string BMbName = stemVec[j];
                shared_ptr<IntFuncTask> request(new IntFuncTask);
                request->asynBind(alignAndScoreMT, pdbA, strandVec[j], pbpl, patl, jid,
                    scoreMap, normedScMap, emptyBMS, emptyBMS, BMaName, BMbName, outPath, outfileAln, outfilePDB);
                jid++;
                thrPool->addTask(request);
            }
            while(true) {
                sleep(1);
                if(thrPool->getTaskCount() == 0) {
                    break;
                }
            }
            cout<<"[Info][Main] Aligment for Motif "<<stemA<<" complete, writing scores."<<endl;
            ofstream scout;
            string scSuf = "_score-" + stemA + ".csv"; 
            string scOutFile = outfileSc.string() + scSuf;
            scout.open(scOutFile, ios::out);
            if (! scout.is_open())
            {
                throw "[Error][Main] Fail to open file " + scOutFile;
            }
            scout << "MotifA,MotifB,score,normedScore" <<endl;
            for(int j=0;j<lv;j++) {
                string BMbName = stemVec[j];
                array<string,2> key{BMaName, BMbName};
                scout << BMaName<<","<<BMbName<<","<<scoreMap->at(key) <<","<<normedScMap->at(key)<<endl;
            }
            scout.close();
            return EXIT_SUCCESS;
        } else {
            throw invalid_argument("pdbB or filelist must be specified when pdbA present");
        }
    }
    if(cmdArgs.specifiedOption("-f")) {
        listFile = cmdArgs.getValue("-f");  // read PDB file names from listfile
        if(cmdArgs.specifiedOption("-d")) {
            cout<<"[Warning][Main] OptionIgnored: Option -d is ignored at the presence of -f"<<endl;
        }
        if(! cmdArgs.specifiedOption("-oa")) {
            outfileAln = "//";
        }
        if(! cmdArgs.specifiedOption("-op")) {
            outfilePDB = "//";
        }
        int nt = 1;
        if(cmdArgs.specifiedOption("-n")) {
            istringstream ss(cmdArgs.getValue("-n"));
            ss >> nt;
        }

        //batch processing
        ifstream input;
        input.open(listFile, ios::in);
        if (! input.is_open())
        {
            throw "[Error][Main] fail to open file " + listFile;
        }
        string s;
        map<int,vector<string> > strandNum2MotifMap;
        shared_ptr<BasePairLib> pbpl(new BasePairLib);
        shared_ptr<AtomLib> patl(new AtomLib);
        while(getline(input,s)) {
            BriqxModule BM0(s, *pbpl, *patl, "BM0", 0, true);
            int nStrand = BM0.getNStrand();
            if(strandNum2MotifMap.contains(nStrand)) {
                strandNum2MotifMap[nStrand].emplace_back(s);
            } else{
                vector<string> vec0{s};
                strandNum2MotifMap.emplace(nStrand, vec0);
            }
            BM0.deepClear();
        }
        input.close();
        cout<<"[Info][Main] List of Motifs Summary:"<<endl;
        cout<<"     nStrand  Population"<<endl;
        for(auto &iter: strandNum2MotifMap) {
            cout<<setw(12)<<iter.first<<setw(12)<<iter.second.size()<<endl;
        }
        shared_ptr<ThreadPool> thrPool(new ThreadPool(nt));
        size_t jid = 0;
        for(auto &iter: strandNum2MotifMap) {
            shared_ptr<map<array<string, 2>, double> > scoreMap(new map<array<string, 2>, double>);
            shared_ptr<map<array<string, 2>, double> > normedScMap(new map<array<string, 2>, double>);
            shared_ptr<BMSelection> emptyBMS(new BMSelection);
            vector<string> stemVec;  // vector of filename w/o suffix
            int lv = iter.second.size();
            for(int i=0;i<lv;i++) {
                filesystem::path path0{iter.second[i]};
                stemVec.emplace_back(path0.stem());
            }
            for(int i=0;i<lv;i++) {  // pre-allocate shared Maps
                string BMaName = stemVec[i];
                for(int j=i+1;j<lv;j++) {
                    string BMbName = stemVec[j];
                    scoreMap->emplace(array<string,2>{BMaName,BMbName}, 0);
                    normedScMap->emplace(array<string,2>{BMaName,BMbName}, 0);
                }
            }
            for(int i=0;i<lv;i++) {
                string BMaName = stemVec[i];
                for(int j=i+1;j<lv;j++) {
                    string BMbName = stemVec[j];
                    shared_ptr<IntFuncTask> request(new IntFuncTask);
                    request->asynBind(alignAndScoreMT, iter.second[i], iter.second[j], pbpl, patl, jid,
                        scoreMap, normedScMap, emptyBMS, emptyBMS, BMaName, BMbName, outPath, outfileAln, outfilePDB);
                    jid++;
                    thrPool->addTask(request);
                }
            }
            while(true) {
                sleep(1);
                if(thrPool->getTaskCount() == 0) {
                    break;
                }
            }
            cout<<"[Info][Main] Aligment for "<<iter.first<<"-strand motifs complete, writing scores."<<endl;
            ofstream scout;
            string scSuf = "_score-" + to_string(iter.first) + "strand.csv"; 
            string scOutFile = outfileSc.string() + scSuf;
            scout.open(scOutFile, ios::out);
            if (! scout.is_open())
            {
                throw "[Error][Main] Fail to open file " + scOutFile;
            }
            scout << "MotifA,MotifB,score,normedScore" <<endl;
            for(int i=0;i<lv;i++) {
                string BMaName = stemVec[i];
                for(int j=i+1;j<lv;j++) {
                    string BMbName = stemVec[j];
                    array<string,2> key{BMaName, BMbName};
                    scout << BMaName<<","<<BMbName<<","<<scoreMap->at(key) <<","<<normedScMap->at(key)<<endl;
                }
            }
            scout.close();
        }
        return EXIT_SUCCESS;
    }
    return EXIT_SUCCESS;
}