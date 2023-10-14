/**
 * @file BMAlign.h
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief  Header file declearing class BMAlign
 * 
 * BMAlign class handles alignment and computes transform paramters of BriqxModules
 * 
 * @version udef
 * @date 2023/08/17
 * 
 * @copyright Copyright (c) 2023 XLAB
 * 
 * modification history :
 * Date:      Version:    Author:
 * Changes:
 */

#ifndef BRIQXMODULE_BMALIGN_H_
#define BRIQXMODULE_BMALIGN_H_
/**
 * @addtogroup BriqxModule
 * @brief BMAlign api description
 * 
 * @{
 */
#include <array>
#include "briqxmodule/BriqxModule.h"
#include "briqxmodule/utils.h"
#include <concepts>
#include <iomanip>

namespace NSPbm
{           
    using namespace NSPgeometry;
    using namespace NSPmodel;
    using namespace std;
    
    /**
     * @brief  BMAlign resolves global alignment bwtween two BriqxModules with an input local alignment called seed
     * 
     */
    template<typename T> // 变量模板, 初始化为对所有类型false
    constexpr bool is_valid_v = false;

    template<> // 变量模板的特化（specialization）, 特化时要把已特化的类型从template参数列表去掉
    constexpr bool is_valid_v<vector<array<RNABase*,2> > > = true;

    template<>
    constexpr bool is_valid_v<vector<array<BasePair*,2> > > = true;

    template<typename T>
    concept IsValid = is_valid_v<T>;
    // template<typename T> requires IsValid<T>  //另一种写法
    template<IsValid T>
    class BMAlign {
        private:
            // bool seedByBasePair;  /** set to true if the instance is seeded by basePair alignment, otherwise set to false */
            // const vector<array<RNABase*, 2> >* baseSeed;  /** local alignment based on RNABase */
            // const vector<array<BasePair*, 2> >* basePairSeed;  /** local alignment based on BasePair */
            const T* seed;
            const BriqxModule* BMa;  /** The first(fixed) BriqxModule object */
            const BriqxModule* BMb;  /** The second(transforming) BriqxModule object */
            vector<array<RNABase*, 2> > alignVec; /** global alignment based on RNAbase */
            vector<array<XYZ,4> > newBList;  /** coordinates of transformed @p BMb bases, represented by FourPseudoAtoms*/
            TransForm tf;  /** Rotational matrix transforms @p BMb to @p newBList */
            XYZ Acog;  /** center-of-geometry of aligned bases of the stational module (module A) */
            XYZ Bcog;  /** center-of-geometry of aligned bases of the transformal module (module B) */

        public:
            BMAlign();
            BMAlign(const BriqxModule& BM1, const BriqxModule& BM2, const T& inSeed); // 类定义内声明成员函数

            const vector<array<RNABase*, 2> >& getAlignVec() const {
                return alignVec;
            }

            const TransForm& getRotTrans() {
                return tf;
            }

            const XYZ getAcog() {
                return Acog;
            }

            const XYZ getBcog() {
                return Bcog;
            }

            int writeAlignment(ofstream& out, double& score, double& normedSc, bool onA = true) const;

            virtual ~BMAlign();
    };

    // template<typename T> requires IsValid<T> // 另一种写法
    template<IsValid T> //  类定义外定义成员函数，template 写法必须和类模板定义的一致
    BMAlign<T>::BMAlign(const BriqxModule& BM1, const BriqxModule& BM2, const T& inSeed) {
	    if(inSeed.size() == 0) return;

        if(BM2.getBaseList().size() == 0) return;

        this->BMa = &BM1;
        this->BMb = &BM2;
        this->seed = &inSeed;
        // this->seedByBasePair = false;

        BriqxModule::resolveTransformByAlign(inSeed, tf, Acog, Bcog);
        BM2.coordTransform(tf, Acog, Bcog, newBList);
        BM1.alignBasesByCoord(BM2, alignVec, newBList);
    }

    // 模板类的实现和声明不能分离，否则连接时会出现undefined reference 错误
    template<IsValid T>
    int BMAlign<T>::writeAlignment(ofstream& out, double& score, double& normedSc, bool onA) const {
        int lav = alignVec.size();
        if(onA) {
            map<RNABase*, RNABase*> a2bMap;
            for(int i=0;i<lav;i++) {
                a2bMap.emplace(get<0>(alignVec[i]), get<1>(alignVec[i]));
            }
            out << "# " << BMa->getBMname() <<"-"<<BMb->getBMname()<<" Base-Base alignment on " << BMa->getBMname() <<endl;
            out << "# Score: " << score <<", Normalized: " << normedSc <<endl;
            out << "# " << setw(6) << "Base1" <<
                           setw(8) << "BaseID" <<
                           setw(8) << "ChainID" <<
                           setw(8) << "Base2" <<
                           setw(8) << "BaseID" <<
                           setw(8) << "ChainID" <<endl;
            auto& bla = BMa->getBaseList();
            int lba = bla.size();
            for(int i=0;i<lba;i++) {
                if(a2bMap.contains(bla[i])) {
                    RNABase* blb = a2bMap[bla[i]];
                    out << setw(8) << bla[i]->baseType <<
                           setw(8) << bla[i]->baseID <<
                           setw(8) << bla[i]->chainID <<
                           setw(8) << blb->baseType <<
                           setw(8) << blb->baseID <<
                           setw(8) << blb->chainID << endl;
                } else {
                    out << setw(8) << bla[i]->baseType <<
                           setw(8) << bla[i]->baseID <<
                           setw(8) << bla[i]->chainID << endl;
                }
            }
        } else {
            map<RNABase*, RNABase*> b2aMap;
            for(int i=0;i<lav;i++) {
                b2aMap.emplace(get<1>(alignVec[i]), get<0>(alignVec[i]));
            }
            out << "# " << BMa->getBMname() <<"-"<<BMb->getBMname()<<" Base-Base alignment on " << BMb->getBMname() <<endl;
            out << "# Score: " << score <<", Normalized: " << normedSc <<endl;
            out << "# " << setw(6) <<"Base1"<<
                           setw(8) <<"BaseID"<<
                           setw(8) <<"ChainID"<<
                           setw(8) <<"Base2"<<
                           setw(8) <<"BaseID"<<
                           setw(8) <<"ChainID"<<endl;
            auto& blb = BMb->getBaseList();
            int lbb = blb.size();
            for(int i=0;i<lbb;i++) {
                if(b2aMap.contains(blb[i])) {
                    RNABase* bla = b2aMap[blb[i]];
                    out << setw(8) << blb[i]->baseType <<
                           setw(8) << blb[i]->baseID <<
                           setw(8) << blb[i]->chainID << 
                           setw(8) << bla->baseType <<
                           setw(8) << bla->baseID <<
                           setw(8) << bla->chainID << endl;
                } else {
                    out << setw(8) << blb[i]->baseType <<
                           setw(8) << blb[i]->baseID <<
                           setw(8) << blb[i]->chainID << endl;
                }
            }
        }
        return EXIT_SUCCESS;
    }

    template<IsValid T>
    BMAlign<T>::~BMAlign() {
        // do nothing
    }

} // namespace NSPbm
#endif /* BRIQXMODULE_BMALIGN_H_ */
/** @} */