/**
 * @file BriqxModule.h
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief  Header file declaring class BriqxModule and BMSelection
 * 
 * BriqxModule is the main class that handles the Briqx modules. BMselection is a supporting class used to locate the Briqx module
 * from a RNAPDB object.
 * 
 * @version udef
 * @date 2023/08/08
 * 
 * @copyright Copyright (c) 2023 XLAB
 * 
 * modification history :
 * Date:      Version:    Author:
 * Changes:
 */

#ifndef BRIQXMODULE_BRIQXMODULE_H_
#define BRIQXMODULE_BRIQXMODULE_H_

/**
 * @addtogroup BriqxModule
 * @brief  BriqxModule api description
 * 
 * Detailed BriqxModule api description
 * 
 * @{
 * 
 */
#include "model/BaseDistanceMatrix.h"
#include "model/StructureModel.h"
#include "geometry/RMSD.h"
#include "model/BasePair.h"
#include "model/BasePairLib.h"
#include <algorithm>
#include <stdexcept>

namespace NSPbm {

    using namespace NSPgeometry;
    using namespace NSPmodel;
    using namespace std;

    /**
     * @brief  The selection class of BriqxModule.
     * 
     * A way to specify a subset of bases from BriqxModule or RNAPDB objects. The selection is based on chains and 
     * bases. Chains are selected alternatively with a @c vector of @c int as 0-based chain indexes (of the
     * original object), or a @c vector of @c string as chain IDs; Bases are selected alternatively with a 2D
     * (axis0->chain, axis1->base) @c vector of @c int as 0-based base indexes of the whole sequence (RNABase baseSeqID), or a 2D
     * @c vector of @c string as base ids in the original PDB (RNABase baseID).
     * 
     */
    class BMSelection {
        private:
            vector<int> chainIndexes;
            vector<vector<int> > resIndexes;
            vector<string> chainIDs;
            vector<vector<string> > chainResIDs;
            bool withChainIndex;  /**< true: use chainIndexes; false: use chainIDs */
            bool withResIndex;  /**< true: use resIndexess; false: use chainResIDs */
        public:
            BMSelection();
            BMSelection(int i);  /**< empty selection, used for NULL input actually */
            BMSelection(const vector<int>& ci, const vector<vector<int> >& cr);
            BMSelection(const vector<int>& ci, const vector<vector<string> >& cr);
            BMSelection(const vector<string>& ci, const vector<vector<int> >& cr);
            BMSelection(const vector<string>& ci, const vector<vector<string> >& cr);
            BMSelection(const string selStr); // TODO: implement selection string parser.

            bool isUseChainIndex() const {
                return withChainIndex;  // 类成员函数内部可以省略this指针，除非函数内定义了同名局部变量
            }
            bool isUseResIndex() const {
                return withResIndex;
            }

            const vector<int>& getChainIndexes() const {  // 返回引用，零次拷贝，但外部应该用普通变量接收（一次拷贝），避免作用域问题
                return chainIndexes;
            }
            const vector<vector<int> >& getResIndexes() const {
                return resIndexes;
            }
            const vector<string>& getChainIDs() const {
                return chainIDs;
            }
            const vector<vector<string> >& getChainResIDs() const {
                return chainResIDs;
            }

            virtual ~BMSelection();
    };

    /**
     * @brief  The main working class of BriqxModule
     * 
     * The class holds the topology of Briqx module by holding the vectors of pointers to the RNAChain and RNABase 
     * objects created by RNAPDB. Functions calculating the distance of DistanceMatrixes, aligning and scoring
     * the alignmments are included.
     * 
     */
    class BriqxModule {
        private:
            string BMname;  /** name of Briqx module */
            vector<RNAChain*> chains;  /** Vector of pointers to RNAchain objects */
                // 指针引用的RNAChain 对象由 new 在堆区创建，不用delete不会自动释放
                // 我们不打算改动 RNAChain 对象，但不保证调用类的操作不会修改 RNAchain 对象
                // 实际上我们也不打算改动指针指向， 但 vector 内置操作要求元素可写，因此不能设为常量指针
            vector<RNABase*> baseList;  /**  Vector of RNABase pointers, size = Nb */
                // 指针引用的RNABase 对象由 new 在堆区创建，不用delete不会自动释放
                // RNABase 对象在计算DM时会被修改
            /**
             * 
             * @brief  Vector of BasePair pointers to BasePairs within the module, size <= Nb(Nb-1)/2 due to filter out
             * of non-contact basepairs.
             */
            vector<BasePair*> basePairList;
            /**
             * 
             * @brief  Vector of BasePair pointers to reversed BasePairs within the module, size=basePairList.size().
             * This list is introduced coz BasePairs are ordered in their DMs.
             */
            vector<BasePair*> revBasePairList; 
            /**
             *
             * @brief Map from 2-elem array of RNAbase* to BasePair*, i.e. from member bases to basepair. Both basepairs 
             * and reversed basepairs are included.
             */
            map<array<RNABase*,2>, BasePair*> bb2bpMap;
            /**
             * 
             * @brief  Map from BasePair to BasePair types in BasePairLib, both BasePair and reversed BasePairs are included.
             */
            map<BasePair*, int> PairTypeMap;
            /**
             * 
             * @brief Map from BasePair to BasePair energies in BasePairLib, both BasePair and reversed BasePairs are included.
             */
            map<BasePair*, double>  PairEneMap;

        public:
            /**
             * @brief  Construct a new empty Briqx Module object
             * 
             */
            BriqxModule();

            /**
             * @brief  Construct a new BriqxModule object from PDB file.
             * 
             * Construct a new Briqx Module object by reading contents of @p pdbFile and get a part of the structure
             * through @p bms . Optionally user can customize the name of the module by @p name . 
             * 
             * @param pdbFile string. Path to the pdb file, required.
             * @param bpl BasePairLib object for basepair type and energy determination, required.
             * @param atl AtomLib object for BasePair determination, required.
             * @param name string. Name of the Briqx module, optional. By default the name of @p pdbFile is used.
             * @param bms BMSelection object. Specifing the range of the module, optional. By default use the whole structure.
             * 
             * @note This constructor allocates new Atom, RNABase, RNAChain and BasePair objects.
             * 
             */
            BriqxModule(const string& pdbFile, BasePairLib& bpl, AtomLib& atl, const string& name = "NewModule", 
                const BMSelection& bms = 0); // BasePair constructor requires non-const Atomlib.

            /**
             * @brief  Construct a new BriqxModule object from given list of chains and list of bases.
             * 
             * @param chains: list of chains; required. 
             * @param baseList: list of bases; required.
             * @param bpl: BasePairLib; required. 
             * @param atl: AtomLib; required. 
             * @param name: name of module; optional, default "NewModule". 
             * @param beLazy: if true, build module in lazy mode, only build @c this->chains and @c this->baseList , else 
             * additionally build basePairList, revBasePairList, PairTypeMap and PaireneMap. Optional, default false
             * @note This constructor may allocate new BasePair objects. It does not check correspondence between chains and
             * baseLists (correctly correlated chains and baseList are expected)
             */
            BriqxModule(const vector<RNAChain*>& chains, const vector<RNABase*>& baseList, BasePairLib& bpl,
                AtomLib& atl, const string& name = "NewModule", bool beLazy = false);
            //TODO: replace the function above with the following two
            BriqxModule(const vector<RNAChain*>& chains, BasePairLib& bpl,
                AtomLib& atl, const string& name = "NewModule", bool beLazy = false);
            BriqxModule(const vector<RNABase*>& baseList, BasePairLib& bpl,
                AtomLib& atl, const string& name = "NewModule", bool beLazy = false);

            /**
             * @brief  Construct a new BriqxModule object by transforming the coordinates of an existing BriqxModule object.
             * 
             * @p tf and @p tv are applied to all atoms within @p bm0 , while BasePair related members are generated only
             * when @p beLazy is false.
             * 
             * @param bm0: The original BriqxModule object. 
             * @param tf: The rotational matrix. 
             * @param tv: The translational vector.
             * @param bpl: BasePairLib; required. 
             * @param atl: AtomLib; required. 
             * @param name: name of module; optional, default use the name of @p bm0. 
             * @param beLazy: if true, build module in lazy mode, only build @c this->chains and @c this->baseList , else 
             * additionally build basePairList, revBasePairList, PairTypeMap and PaireneMap. Optional, default false
             * 
             * @note This constructor allocates new Atom, RNABase, RNAChain and BasePair objects, which should be destructed
             * before destructing objects constructed with this constructor. BriqxModule::deepClear() can be employed for the
             * job.  
             */
            BriqxModule(const BriqxModule& bm0, const TransForm& tf, const XYZ& tv, BasePairLib& bpl,
                AtomLib& atl, const string& name = NULL, bool beLazy = false);

            /**
             * @brief  Set name of the module
             * 
             * @param name string. Name of the module 
             */
            void setBMname(const string& name) {
                BMname = name;
            }

            /**
             * @brief  Get name of the module
             * 
             * @return const sting& BMname the name of the module
             * @note the reference is const and not editable unless BMname is received by a usual string object. 
             */
            const string& getBMname() const {
                return BMname;
            }

            /**
             * @brief  Get ID of the chains in the module 
             * 
             * @return const vector<string>& Reference to a vector of the chain IDs
             * @note   The referenced vector from getChainIDs is allocated by @c new operator,
             *         thus should be deleted manually.
             */
            const vector<string>& getChainIDs() const {
                vector<string>* chainIDs = new vector<string>;
                int len = chains.size();
                for(int i = 0; i < len; i++) { // oldest ways result fastest iteration
                    chainIDs->emplace_back(chains[i]->getChainID());
                }
                return *chainIDs;
            }   

            const vector<RNABase*>& getBaseList() const {
                return this->baseList;
            }

            /**
             * @brief  Get the BasePairList within the BriqxModule object
             * 
             * @return const vector<BasePair*>& 
             */
            const vector<BasePair*>& getBasePairList() const {
                return this->basePairList;
            }

            /**
             * @brief  Get the revBasePairList within the BriqxModule object
             * 
             * @return const vector<BasePair*>& 
             */
            const vector<BasePair*>& getRevBasePairList() const {
                return this->revBasePairList;
            }

            const map<array<RNABase*,2>, BasePair*>& getBb2bpMap () const {
                return bb2bpMap;
            }

            const map<BasePair*, double>& getBasePairEneMap() const {
                return PairEneMap;
            }

            const map<BasePair*, int>& getBasePairTypeMap() const {
                return PairTypeMap;
            }

            /**
             * @brief  Get a list of BasePairs within the module according to a given @p bms
             * 
             * @param bms A BMSelection object specifiing the range of bases to take account of, optional.
             *            By default return all BasePairs within the module.
             * @param withReversed A flag specifying whether return reversed BasePairList, if true, the basepairs from
             * revBasePairList is included, otherwise only basepairs from basePairList are included.
             * @return vector<BasePair*>&  A reference to a vector of BasePairs.
             * @note  The returned vector is created by New operator and should be released manually.
             */
            vector<BasePair*>& getBasePairList(const BMSelection& bms, const bool withReversed = false) const;

            /**
             * @brief  Calculate a matrix of Distances between base Distance Matrixes (DDM) to another module.
             * 
             * For each basepair in current module, calculate the Distance bwteen its base Distance Matrix and that of
             * each basepair from another module, which composes the Matrix of DDM, with row axis being basepairs (
             * reverse basepairs included) from the current module and the column axis being basepairs (reverse
             * basepairs ignored) from the other module.
             * 
             * @param [in] other: Reference to a BriqxModule object.
             *      The other module.
             * @param [out] DDMMatrix: Matrix of the calculated DDMs, initialized as a map mapping a 2-elem array
             *      denoting BasePairs from @c current and @p other module respectively, to the DDM value.
             * 
             * @return int, code of successfulness.
             */
            int calcDDMMatrix(const BriqxModule& other, map<array<BasePair*, 2>, double>& DDMMatrix) const;

            /**
             * @brief  Sort given map by value.
             * 
             * Sort given map by value, return a vector composed of sorted key-value pairs. 
             * 
             * @param [in] unsortMap:
             * @param [in] order: order of sorting, accepts "desc", "asc", default "desc" corresponds to descending order
             * @param [out] sortedVec: 
             * 
             * @return int: code of successfulness.
             * @note replaced by utils::sortMapByValue
             */
            static int sortMapByValue(const map<array<BasePair*, 2>, double>& unsortMap,
                vector<pair<array<BasePair*, 2>,double> >& sortedVec, string order = "desc");

            /**
             * @brief  Resolve the transform to a BriqxModule according to aligned BasePairs ( @p alignVec ).
             * 
             * The transform includes a rotational matrix @p tf and a translational vector @p tv .
             * 
             * @param [in] alignVec: vector of pairs of aligned BasePairs. In each element, the first BasePair is from the first 
             * module and the second BasePair is from the second one. 
             * @param [out] tf: A TransForm object that holds rotational matrix, which when applied will rotate the second module
             * to the first one.
             * @param [out] tv: An XYZ object that holds the translation vector, which when applied will translate the second
             * module to the first one.
             * @return int : success code.
             */
            static int resolveTransformByAlign(const vector<array<BasePair*,2> >& alignVec, TransForm& tf, XYZ& tv);

            /**
             * @brief  Resolve the transform to a BriqxModule according to aligned Bases ( @p alignVec ).
             * 
             * The transform includes a rotational matrix @p tf and a translational vector @p tv .
             * 
             * @param [in] alignVec: vector of pairs of aligned RNABases. In each element, the first RNABase is from the first 
             * module and the second RNABase is from the second one. 
             * @param [out] tf: A TransForm object that holds rotational matrix, which when applied will rotate the second module
             * to the first one.
             * @param [out] tv: An XYZ object that holds the translation vector, which when applied will translate the second
             * module to the first one.
             * @return int : success code.
             * @note Function called by BMAlign, for usual taskes invoked through BMAlign is recommended.
             */
            static int resolveTransformByAlign(const vector<array<RNABase*,2> >& alignVec, TransForm& tf, XYZ& tv);

            /**
             * @brief  Apply transform (TransForm @p tm  and translational vector @p tv ) to @c this BriqxModule, return
             * transformed FourPseudoAtoms coordinates of bases in @c this BriqxModule. 
             * 
             * @param [in] tf: A TransForm object that rotate @c this module to the orientation of another module.
             * @param [in] tv: A vector that translates @c this module to another module (by center of geometry).
             * @param [out] tfBasePos: vector of 4-element arrays corresponding to the four pseudoAtoms in each base.
             * Each element of the arrays stores the transformed coordinates of the pesudoAtoms.
             * @return int: success code.
             * @note Function called by BMAlign, for usual taskes invoked through BMAlign is recommended.
             */
            int coordTransform(const TransForm& tf, const XYZ& tv, vector<array<XYZ,4> >& tfBasePos) const;

            /**
             * @brief  Resolve base-base alignment to another BriqxModule according to base distance. The base distance
             * is defined as the mean of distances between corresponding FourPseudoAtoms. Bases are aligned to the nearest
             * base in @p other module.
             * 
             * @param [in] other: The other BriqxModule. 
             * @param [in] dCutoff: Cutoff of base-base distance to treat two bases as aligned.
             * @param [in] BasePos: New position of the bases of @p other . The alignment evaluation will base on
             * it if it is not empty. Optional, default empty.
             * @param [out] alignVec: <array<RNABase*,2> > vector of pairs of aligned bases.
             * @return int: success code.
             * @note Function called by BMAlign, for usual taskes invoked through BMAlign is recommended.
             */
            int alignBasesByCoord(const BriqxModule& other, vector<array<RNABase*, 2> >& alignVec,
            const vector<array<XYZ,4> >& BasePos = vector<array<XYZ,4> >{}, const double& distMax = 1.0) const;

            /**
             * @brief  Print object in PDB format to @p out.
             * 
             * @param out: Output stream.
             */
            void printPDBFormat(ofstream& out) const;

            /**
             * @brief  Delete all allocated memory that referenced by pointers within the BriqxModule object.
             * RNABase, BasePair and RNAChain objects will be deleted. This function is expected to be called
             * only before destruction. 
             * 
             * @note Class members holding the pointers to these objects are not removed (they are expected
             * to be deleted during destruction)
             * @warning Call with caution that the cleared objects are no longer used by any other instances.
             * 
             */
            void deepClear();

            /**
             * @brief  Delete BasePair objects referenced by basePairList and revBasePairList. This function is
             * expected to be called only before destruction.
             * 
             * @note Class members holding the pointers to these objects are not removed (they are expected
             * to be deleted during destruction)
             * @warning Call with caution that the cleared objects are no longer used by any other instances.
             * 
             */
            void clearBasePairs();
            
            virtual ~BriqxModule();
    };
} // NSPbm

#endif /* BRIQXMODULE_BRIQXMODULE_H_ */
/** @} */