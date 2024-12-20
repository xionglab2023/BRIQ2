/**
 * @file BMScore.h
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief  Header file for class BMScore
 *
 * BMScore handles scoring of an BriqxModule alignment
 *
 * @version udef
 * @date 2023/08/19
 *
 * @copyright Copyright (c) 2023 XLAB
 *
 * modification history :
 * Date:      Version:    Author:
 * Changes:
 */

#ifndef BRIQXMODULE_BMSCORE_H_
#define BRIQXMODULE_BMSCORE_H_

/**
 * @addtogroup BriqxModule
 * @brief BMScore api description
 *
 * @{
*/
#include "briqxmodule/BriqxModule.h"
#include "briqxmodule/BMAlign.h"

namespace NSPbm
{
    using namespace NSPgeometry;
    using namespace NSPmodel;
    using namespace std;

    /**
     * @brief  BMScore objects holds the weight and DDM info for scoring, and scores a given alignment or a list of pair of BasePairs
     *
     */
    class BMScore {
        private:
            const BriqxModule* BMa;
            const BriqxModule* BMb;
            map<BasePair*, double> bp2wMapa;  /** Map from @c BMa BasePair* to its score weight */ //weight formula is defined in BMscore
            map<BasePair*, double> bp2wMapb;  /** Map from @c BMb BasePair* to its score weight */
            const map<array<BasePair*, 2>, double>* DDMMap;  /** Map from pair of basepairs (one bp from @c BMa and the other from @c BMb) to DDM */
                /* DDM is determined by @c BMa and @c BMb, not necessary to get one for each BMscore instance, but also not suitable
                * to make it static coz the class may handle different BMa and BMb */

            /**
             * @brief Unify input alignVec to vector of pairs (2-elem array) of BasePairs.
             *
             * If receive a vector of pairs of BasePairs, return a deep copy of the vector, else if receive a vector of pairs of RNABases, explore
             * possible BasePairs and return the BasePairs in pairs.
             *
             * @param inAlignVec: Input align vector
             * @return vector<array<BasePair*,2> >&
             * @note Called by score(). The returned vector need to be relased manually.
             */
            const vector<array<BasePair*,2> >* unifyAlignVec(const vector<array<BasePair*,2> >& inAlignVec);
            const vector<array<BasePair*,2> >* unifyAlignVec(const vector<array<RNABase*,2> >& inAlignVec);
        public:
            BMScore();
            BMScore(const BriqxModule& BM1, const BriqxModule& BM2, const map<array<BasePair*, 2>, double>& DDMMatrix);
            /**
             * @brief  Score a given alignment. Function template, instancialized with vector of pair of bases or basepairs
             *
             * @tparam T:  vector, contains pairs of RNABase*s or BasePair*s
             * @param inAlignVec: Input align vector in type T
             * @return double: The score of the given alignment.
             */
            template<typename T> requires IsValid<T>
            double score(const T& inAlignVec);

            /**
             * @brief  Calculate the score of ideal alignment (thus called ideal score)
             *
             * An ideal alignment is defined as a full match between BMa and BMb, i.e. every base(basepair) matches, thus
             * all weights are counted and all DDM=0.
             *
             * @return double
             */
            double idealScore(bool bmf = false);
            double idealScore2(); // scoring according to strands restricted to the length of the shorter strand

            /**
             * @brief  Sort BasePair pairs by score.
             *
             * @param [in] bppList: Vector of pairs of BasePairs
             * @param [out] sortedScore:  Vector of pairs of BasePairs that sorted by score.
             * @param [in] desc: If true, sort in descending order, else ascending.
             * @return int: success code
             */
            int sortBppByScore(const vector<array<BasePair*,2> >& bppList, vector<pair<array<BasePair*,2>, double> >& sortedScore,
                bool desc=true);
            int sortBppByScore(vector<pair<array<BasePair*,2>, double> >& sortedScore, bool desc=true);  // sort BasePairs in DDMMap

            /**
             * @brief  Transfer given energy value to weight.
             *
             * @param ene: BasePair energy.
             * @return double: The weight in bpp score.
             */
            static double getWeight(double ene) {
                // return 1/(1+exp(ene/2+5)); // if ene in kcal/mol
                // return 1/(1+exp(ene*3+6)); // if ene in BRIQX unit, 中心于-2 BRIQX UNIT，饱和于-4 BRIQX UNIT
                return 1/(1+exp(ene*2.4+3)); // if ene in BRIQX unit, 中心于-1.25 BRIQX UNIT，饱和于-2.5 BRIQX UNIT
            }

            /**
             * @brief  Calculate the socre of a single pair of BasePair (bpp).
             *
             * @param wa: Weight of BasePair a in the pair.
             * @param wb: Weight of BasePair b in the pair.
             * @param ddm: DDM between BasePair a and b.
             * @return double: The score of the BasePair pair.
             */
            static double bppScore(double wa, double wb, double ddm) {
                // return (wa+wb)/(ddm+1)*0.5;  // round 5-7
                // return max((wa+wb)*0.5*(1-ddm/2),0.0);  // round 8
                // return max((wa+wb)*0.5*(1-ddm/3),0.0);  // round 9
                // return max((wa+wb)*0.5*(1-ddm*ddm/4),0.0);  //round 10
                // return (wa+wb)*0.5*(exp(-ddm*ddm/0.98));  // round 11
                return (wa+wb)*0.5*(exp(-ddm*ddm/2));  // round 12, current
                // return (wa+wb)*0.5*(exp(-ddm*ddm/2.88));  // round 13
                // return (wa+wb)*0.5*(exp(-ddm*ddm/3.92));  // round 14
                // return (wa+wb)*0.5*(exp(-ddm*ddm/1.62));  // round 15
            }
            
            #ifdef DEBUG
            /**
             * @brief Print bp2wMapa and bp2wMapb to @param out in human-friendly form
             * 
             * @param out: Output stream
            */
            void printBasePairWeight(ostream& out);
            #endif

            virtual ~BMScore();
    };

    template<typename T> requires IsValid<T>
    double BMScore::score(const T& inAlignVec) {
        const vector<array<BasePair*,2> >* alignVec = unifyAlignVec(inAlignVec);
        int lav = alignVec->size();
        double tot = 0;
        for(int i=0;i<lav;i++) {
            double curBppSc = bppScore(
                bp2wMapa.at((*alignVec)[i][0]), bp2wMapb.at((*alignVec)[i][1]), DDMMap->at((*alignVec)[i]));
            #ifdef DEBUG
            string outstr = "[DEBUG][Info][score] weight " + (*alignVec)[i][0]->print() + ": " + \
                to_string(bp2wMapa.at((*alignVec)[i][0])) + "\n";
            cout << outstr;
            outstr = "[DEBUG][Info][score] weight " + (*alignVec)[i][1]->print() + ": " + \
                to_string(bp2wMapb.at((*alignVec)[i][1])) + "\n";
            cout << outstr;
            outstr = "[DEBUG][Info][score] DDM: " + to_string(DDMMap->at((*alignVec)[i])) + "\n";
            cout << outstr;
            outstr = "[DEBUG][Info][score] bppScore: " + to_string(curBppSc) + "\n";
            cout << outstr;
            #endif
            tot += curBppSc;
        }
        delete alignVec;
        return tot;
    }


} // namespace NSPbm

#endif /* BRIQXMODULE_BMSCORE_H_ */
/** @} */
