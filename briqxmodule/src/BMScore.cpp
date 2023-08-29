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
        int lbpa = bpa.size();
        for(int i=0;i<lbpa;i++) {
            bp2wMapa.emplace(bpa[i], getWeight(bpEnea.at(bpa[i])));
            bp2wMapa.emplace(rbpa[i], getWeight(bpEnea.at(rbpa[i])));
        }
        const vector<BasePair*>& bpb = BM2.getBasePairList();
        const vector<BasePair*>& rbpb = BM2.getRevBasePairList();
        const map<BasePair*, double>& bpEneb = BM2.getBasePairEneMap();
        int lbpb = bpb.size();
        for(int i=0;i<lbpb;i++) {
            bp2wMapb.emplace(bpb[i], getWeight(bpEneb.at(bpb[i])));
            bp2wMapb.emplace(rbpb[i], getWeight(bpEneb.at(rbpb[i])));
        }
    }

    const vector<array<BasePair*,2> >& BMScore::unifyAlignVec(const vector<array<BasePair*, 2> >& inAlignVec) {
        auto* out = new vector<array<BasePair*,2> >{inAlignVec};
        return *out;
    }

    const vector<array<BasePair*,2> >& BMScore::unifyAlignVec(const vector<array<RNABase*,2> >& inAlignVec) {
        auto* out = new vector<array<BasePair*,2> >{};
        auto& bb2bpMapA = BMa->getBb2bpMap();
        auto& bb2bpMapB = BMb->getBb2bpMap();
        int lav = inAlignVec.size();
        for(int i=0;i<lav;i++) {
            RNABase* bA1 = inAlignVec[i][0];
            RNABase* bB1 = inAlignVec[i][1];
            for(int j=i+1;j<lav;j++) {
                RNABase* bA2 = inAlignVec[j][0];
                RNABase* bB2 = inAlignVec[j][1];
                array<RNABase*, 2> arrBbA{bA1, bA2};  // array of base-base, the alignment is ordered, no reverse to handle
                array<RNABase*, 2> arrBbB{bB1, bB2};
                if(bb2bpMapB.contains(arrBbB) && bb2bpMapA.contains(arrBbA)){
                    out->emplace_back(array<BasePair*,2>{bb2bpMapA.at(arrBbA), bb2bpMapB.at(arrBbB)});
                } else {
                    cout << "Warning: " << BMa->getBMname() << ":" << bA1->baseType << bA1->baseID << "-" << \
                    bA2->baseType<<bA2->baseID << " or " \
                    << BMb->getBMname() << ":" << bB1->baseType << bB1->baseID << "-" \
                    << bB2->baseType << bB2->baseID << " BasePair not qualified" << endl;
                }
            }
        }
        return *out;
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
        for (auto& iter : bp2wMapa) {
            score += iter.second;
        }
        for (auto& iter : bp2wMapb) {
            score += iter.second;
        }
        return score/2;
    }

    BMScore::~BMScore(){
        
    }

} // namespace NSPbm
