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

    using namespace NSPtools;
    using namespace NSPbm;
    using namespace std;

    /**
     * @brief Workhorse function to align and calculate alignment score for two BriqxModules 
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
        const string& outfileAln = "align.aln", const string& outfilePDB = "alignedB.pdb"){
        
        const double DDMCutoff = 0.7;
        const double iterStopAt = 1e-6;

        BriqxModule BMa(BMaPDB, bpl, atl, BMaName, BMaSel);
        BriqxModule BMb(BMbPDB, bpl, atl, BMbName, BMbSel);

        map<array<BasePair*, 2>, double> DDMMatrix{};

        BMa.calcDDMMatrix(BMb, DDMMatrix);

        BMScore scer(BMa, BMb, DDMMatrix);

        vector<pair<array<BasePair*, 2>, double> > sortedDDM{};

        scer.sortBppByScore(sortedDDM);

        int i = 0;
        auto* p1 = &sortedDDM[i];
        double bestScore = 0;
        BMAlign<vector<array<RNABase*,2> > >* bestAlign;
        while(p1->second <= DDMCutoff) {
            BMAlign<vector<array<BasePair*,2> > > iniAlign(BMa, BMb, vector<array<BasePair*,2> >{p1->first});
            BMAlign<vector<array<RNABase*,2> > >* curAlign = NULL, *lastAlign, *nextAlign;
            double curSc = 0, lastSc;
            double nextSc = scer.score(iniAlign.getAlignVec());
            if(curSc-nextSc < iterStopAt) {
                lastSc = curSc;
                curSc = nextSc;
                nextAlign = new BMAlign<vector<array<RNABase*,2> > >(BMa, BMb, iniAlign.getAlignVec());
                nextSc = scer.score(nextAlign->getAlignVec());
            }
            while(curSc-nextSc < iterStopAt) {
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
                bestAlign = curAlign;
            } else {
                delete curAlign;
            } 
            try {   
                p1 = &sortedDDM.at(++i);
            } catch(out_of_range) {
                break;
            }
        }
        
        score = bestScore;
        alignVec = bestAlign->getAlignVec();
        double idealSc = scer.idealScore();
        normedSc = score/idealSc;

        cout<< BMa.getBMname() <<"-"<<BMb.getBMname()<<" Base-Base alignment completed" <<endl;
        cout<<"Alignment score:" << to_string(score) << endl;
        cout<<"Normed score:" << to_string(normedSc) << endl;

        if(outPath != "") {
            cout<<"Write alignment on "<<BMa.getBMname()<<" to " <<outPath<<"/"<<outfileAln<<"_onA"<<endl;
            ofstream output;
            string filePath = outPath+"/"+outfileAln+"_onA";
            output.open(filePath, ios::out);
            bestAlign->writeAlignment(output, score, normedSc, true);
            output.close();
            cout<<"Write alignment on "<<BMb.getBMname()<<" to " <<outPath<<"/"<<outfileAln<<"_onB"<<endl;
            filePath = outPath+"/"+outfileAln+"_onB";
            output.open(filePath, ios::out);
            bestAlign->writeAlignment(output, score, normedSc, false);
            output.close();
            cout<<"Write transformed "<<BMb.getBMname()<<" to " <<outPath<<"/"<<outfilePDB<<endl;
            output.open(filePath, ios::out);
            auto BMnew = BriqxModule(BMb, bestAlign->getRotTrans(), bestAlign->getTranslateVec(), bpl, atl, NULL, true);
            BMnew.printPDBFormat(output);
            BMnew.deepClear();
        }

        BMa.deepClear();
        BMb.deepClear();
        delete bestAlign;
        return EXIT_SUCCESS;
    }

    void printHelp() {
        cout << "Usage:" <<endl;
        cout << "bmalign -a PDBA -b PDBB -o outPath -op outPDBFileName -oa outAlignmentFileName" << endl;
    }

    /**
     * @brief  Main function of the Program.
     * 
     * @param Options accepted:
     * @param -h: print help
     * @param -a: path to module A PDB file
     * @param -b: Path to module B PDB file 
     * @param -o: Output path
     * @param -op: filename (no path included) of output PDB file recording aligned module B
     * @param -oa: filename (no path included) of output alignment files (based on module A and module B, respectively).
     * @param Options not implemented:
     * @param -f: Read input PDB list from given file.
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
        string listfile, pdbA, pdbB, sel;
        string outPath{};
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
        if(cmdArgs.specifiedOption("-sel")) {
            sel = cmdArgs.getValue("-sel");   // parser for sel string has not been implemented.
        }
        if(cmdArgs.specifiedOption("-f")) {
            listfile = cmdArgs.getValue("-f");  // read PDB file names from listfile
            //batch processing
        } else if(cmdArgs.specifiedOption("-a")) {
            pdbA = cmdArgs.getValue("-a");  // -a and -b must be specified at the same time
            if(cmdArgs.specifiedOption("-b")) {
                pdbB = cmdArgs.getValue("-b");
            } else {
                throw invalid_argument("pdbB must be specified when pdbA present");
            }
            BasePairLib bpl;
            AtomLib atl;
            alignAndScore(pdbA, pdbB, bpl, atl, score, alignVec, normedSc, 0, 0,  // for now use NULL(i.e. 0) BMSelection
                "BMa", "BMb", outPath, outfileAln, outfilePDB);
        }
        return EXIT_SUCCESS;
    }