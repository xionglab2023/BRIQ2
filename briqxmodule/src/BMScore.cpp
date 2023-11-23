/**
 * @file BMScore.cpp
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief
 * @version udef
 * @date 2023/08/19
 *
 * @copyright Copyright (c) 2023 XLAB
 *
 * modification history :
 * Date:      Version:    Author:
 * Changes:
 */

#include "briqxmodule/BMScore.h"
#include "model/StructureModel.h"
#include <array>

namespace NSPbm
{
    using namespace NSPgeometry;
    using namespace NSPmodel;
    using namespace std;

    BMScore::BMScore() {
        BMa = NULL;
        BMb = NULL;
        DDMMap = NULL;
    }

    BMScore::BMScore(const BriqxModule& BM1, const BriqxModule& BM2, const map<array<BasePair*, 2>, double>& DDMMatrix)
    {
        BMa = &BM1;
        BMb = &BM2;
        DDMMap = &DDMMatrix;
        const vector<BasePair*>& bpa = BM1.getBasePairList();
        const vector<BasePair*>& rbpa = BM1.getRevBasePairList();
        const map<BasePair*, double>& bpEnea = BM1.getBasePairEneMap();
        // reweight terminal BasePairs
        vector<RNABase*> termBaseVec;
        const auto& stda = BM1.getStrands();
        int lsa = stda.size();
        for(int i=0; i<lsa; i++) {
            const auto curBaseList = stda[i]->getBaseList();
            int stdLen = stda[i]->getChainLength();
            if(stdLen > 6) {
                termBaseVec.emplace_back(curBaseList[0]);
                termBaseVec.emplace_back(curBaseList[1]);
                termBaseVec.emplace_back(curBaseList[stdLen-1]);
                termBaseVec.emplace_back(curBaseList[stdLen-2]);
            } else if(stdLen > 4) {
                termBaseVec.emplace_back(curBaseList[0]);
                termBaseVec.emplace_back(curBaseList[stdLen-1]);
            }
        }
        set<BasePair*> termBpSet;
        int ltbv = termBaseVec.size();
        const auto& bb2bpMapa = BM1.getBb2bpMap();
        for(int i=0; i<ltbv; i++) {
            for(int j=i+1; j<ltbv; j++) {
                array<RNABase*, 2> key = array<RNABase*, 2>{termBaseVec[i], termBaseVec[j]};
                array<RNABase*, 2> rkey = array<RNABase*, 2>{termBaseVec[j], termBaseVec[i]};
                if(bb2bpMapa.count(key) != 0) {
                    termBpSet.insert(bb2bpMapa.at(key));
                    termBpSet.insert(bb2bpMapa.at(rkey));
                }
            }
        }
        //reweight terminal basepairs done
        double termFactor = 1;
        int lbpa = bpa.size();
        for(int i=0;i<lbpa;i++) {
            if(termBpSet.count(bpa[i]) != 0) {
                termFactor = 0.5;
            } else {
                termFactor = 1;
            }
            bp2wMapa.emplace(bpa[i], termFactor*getWeight(bpEnea.at(bpa[i])));
            bp2wMapa.emplace(rbpa[i], termFactor*getWeight(bpEnea.at(rbpa[i])));
        }

        const vector<BasePair*>& bpb = BM2.getBasePairList();
        const vector<BasePair*>& rbpb = BM2.getRevBasePairList();
        const map<BasePair*, double>& bpEneb = BM2.getBasePairEneMap();
        // reweight terminal BasePairs
        termBaseVec.clear();
        const auto& stdb = BM2.getStrands();
        int lsb = stdb.size();
        for(int i=0; i<lsb; i++) {
            const auto curBaseList = stdb[i]->getBaseList();
            int stdLen = stdb[i]->getChainLength();
            if(stdLen > 6) {
                termBaseVec.emplace_back(curBaseList[0]);
                termBaseVec.emplace_back(curBaseList[1]);
                termBaseVec.emplace_back(curBaseList[stdLen-1]);
                termBaseVec.emplace_back(curBaseList[stdLen-2]);
            } else if(stdLen > 4) {
                termBaseVec.emplace_back(curBaseList[0]);
                termBaseVec.emplace_back(curBaseList[stdLen-1]);
            }
        }
        termBpSet.clear();
        ltbv = termBaseVec.size();
        const auto& bb2bpMapb = BM2.getBb2bpMap();
        for(int i=0; i<ltbv; i++) {
            for(int j=i+1; j<ltbv; j++) {
                array<RNABase*, 2> key = array<RNABase*, 2>{termBaseVec[i], termBaseVec[j]};
                array<RNABase*, 2> rkey = array<RNABase*, 2>{termBaseVec[j], termBaseVec[i]};
                if(bb2bpMapb.count(key) != 0) {
                    termBpSet.insert(bb2bpMapb.at(key));
                    termBpSet.insert(bb2bpMapb.at(rkey));
                }
            }
        }
        //reweight terminal basepairs done
        int lbpb = bpb.size();
        for(int i=0;i<lbpb;i++) {
            if(termBpSet.count(bpb[i]) != 0) {
                termFactor = 0.5;
            } else {
                termFactor = 1;
            }
            bp2wMapb.emplace(bpb[i], termFactor*getWeight(bpEneb.at(bpb[i])));
            bp2wMapb.emplace(rbpb[i], termFactor*getWeight(bpEneb.at(rbpb[i])));
        }
    }

    const vector<array<BasePair*,2> >* BMScore::unifyAlignVec(const vector<array<BasePair*, 2> >& inAlignVec) {
        auto* out = new vector<array<BasePair*,2> >{inAlignVec};
        return out;
    }

    const vector<array<BasePair*,2> >* BMScore::unifyAlignVec(const vector<array<RNABase*,2> >& inAlignVec) {
        auto* out = new vector<array<BasePair*,2> >{};
        auto& bb2bpMapA = BMa->getBb2bpMap();
        auto& bb2bpMapB = BMb->getBb2bpMap();
        auto& blB = BMb->getBaseList();
        int lav = inAlignVec.size();
        for(int i=0;i<lav;i++) {
            RNABase* bA1 = inAlignVec[i][0];
            RNABase* bB1 = inAlignVec[i][1];
            for(int j=i+1;j<lav;j++) {
                RNABase* bA2 = inAlignVec[j][0];
                RNABase* bB2 = inAlignVec[j][1];
                /* Ensure base order in basepair:
                 *   For module A, both forward and reverse basepairs are accepted;
                 *   For module B, only accrpt forward basepairs
                */
                array<RNABase*, 2> arrBbA, arrBbB;
                if(find(blB.begin(), blB.end(), bB1) - find(blB.begin(), blB.end(), bB2) < 0) {
                    // bB1-bB2 is forward basepair
                    arrBbA[0] = bA1;
                    arrBbA[1] = bA2;
                    arrBbB[0] = bB1;
                    arrBbB[1] = bB2;
                } else {
                    // bB2-bB1 is forward basepair
                    arrBbA[0] = bA2;
                    arrBbA[1] = bA1;
                    arrBbB[0] = bB2;
                    arrBbB[1] = bB1;
                }
                if(bb2bpMapB.contains(arrBbB) && bb2bpMapA.contains(arrBbA)){
                    out->emplace_back(array<BasePair*,2>{bb2bpMapA.at(arrBbA), bb2bpMapB.at(arrBbB)});
                // } else {
                    // #ifdef DEBUG
                    // string outstr = "[DEBUG][Warning] " + BMa->getBMname() + ":" + bA1->baseType + bA1->baseID + "-" \
                                    //  + bA2->baseType + bA2->baseID + " or " + \
                                    //  BMb->getBMname() + ":" + bB1->baseType + bB1->baseID + "-" \
                                    //  + bB2->baseType + bB2->baseID + " BasePair not qualified" + "\n";
                    // cout << outstr;
                    // #endif
                }
            }
        }
        return out;
    }

    int BMScore::sortBppByScore(const vector<array<BasePair*,2> >& bppList, vector<pair<array<BasePair*,2>, double> >& sortedScore,
        bool desc) {
        int ll = bppList.size();
        map<array<BasePair*,2>, double> tmpMap;
        for(int i=0;i<ll;i++) {
            double score = bppScore(bp2wMapa.at(bppList[i][0]), bp2wMapb.at(bppList[i][1]), DDMMap->at(bppList[i]));
            tmpMap.emplace(bppList[i], score);
        }
        return utils::sortMapByValue<array<BasePair*,2>, double>(tmpMap, sortedScore, desc);
    }

    int BMScore::sortBppByScore(vector<pair<array<BasePair*,2>, double> >& sortedScore, bool desc) {
        auto tmpMap{*DDMMap};
        for(auto& iter : tmpMap) {
            double score = bppScore(bp2wMapa.at(iter.first[0]), bp2wMapb.at(iter.first[1]), iter.second);
            iter.second = score;
        }
        return utils::sortMapByValue<array<BasePair*,2>, double>(tmpMap, sortedScore, desc);
    }

    double BMScore::idealScore() {
        double score = 0;
        int lA = BMa->getBaseList().size();
        int lB = BMb->getBaseList().size();
        if(lA == lB || lA -lB > lB* LDIFF_RATIO_TOLERATE || lB - lA > lA*LDIFF_RATIO_TOLERATE) {
            // Ideally all matches are between forward basepairs
            auto& bplA = BMa->getBasePairList();
            int lbA = bplA.size();
            for(int i=0;i<lbA;i++) {
                score += bp2wMapa.at(bplA[i]);
            }
            auto& bplB = BMb->getBasePairList();
            int lbB = bplB.size();
            for(int i=0;i<lbB;i++) {
                score += bp2wMapb.at(bplB[i]);
            }
            // for (auto& iter : bp2wMapa) {
            //     score += iter.second;
            // }
            // for (auto& iter : bp2wMapb) {
            //     score += iter.second;
            // }
            return score/2;
        } else {
            if(lA > lB) {
                auto& bplB = BMb->getBasePairList();
                int lbB = bplB.size();
                for(int i=0;i<lbB;i++) {
                    score += bp2wMapb.at(bplB[i]);
                }
            } else {
                auto& bplA = BMa->getBasePairList();
                int lbA = bplA.size();
                for(int i=0;i<lbA;i++) {
                    score += bp2wMapa.at(bplA[i]);
                }
            }
            return score;
        }
    }

    double BMScore::idealScore2() {
        double score=0;
        return score;
    }

    #ifdef DEBUG
    void BMScore::printBasePairWeight(ostream& out) {
        out<<"[DEBUG][Info] Weight for BasePairs of "<<BMa->getBMname()<<":"<<endl;
        auto& bpList = BMa->getBasePairList();
        auto& rbpList = BMa->getRevBasePairList();
        int lbp = bpList.size();
        for(int i=0; i<lbp;i++) {
            auto* curBP = bpList[i];
            auto* currBP = rbpList[i];
            out << curBP->print() << ":  " << bp2wMapa.at(curBP) <<endl;
            out << currBP->print() << ":  " << bp2wMapa.at(currBP) <<endl;
        }
        out<<"[DEBUG][Info] Weight for BasePairs of "<<BMb->getBMname()<<":"<<endl;
        auto& bpList2 = BMb->getBasePairList();
        auto& rbpList2 = BMb->getRevBasePairList();
        lbp = bpList2.size();
        for(int i=0; i<lbp;i++) {
            auto* curBP = bpList2[i];
            auto* currBP = rbpList2[i];
            out << curBP->print() << ":  " << bp2wMapb.at(curBP) <<endl;
            out << currBP->print() << ":  " << bp2wMapb.at(currBP) <<endl;
        }
    }
    #endif

    BMScore::~BMScore(){

    }

} // namespace NSPbm
