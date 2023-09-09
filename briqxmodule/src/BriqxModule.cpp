/**
 * @file BriqxModule.cpp
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief  Implementation of BMSelection and BriqxModule classes
 * @version udef
 * @date: 2023/08/08
 *
 * @copyright Copyright (c) 2023 XLAB
 *
 * modification history :
 * Date:      Version:    Author:
 * Changes:
 */

#include "briqxmodule/BriqxModule.h"
#include "briqxmodule/utils.h"
#include <array>

namespace NSPbm {

BMSelection::BMSelection() {
    this->chainIndexes.clear();
    this->chainIDs.clear();
}

BMSelection::BMSelection(int i) {
    this->chainIndexes.clear();
    this->chainIDs.clear();
}

BMSelection::BMSelection(const vector<int>& ci, const vector<vector<int> >& cr) {
    if(ci.size() != cr.size()) {
        throw invalid_argument("Inconssitent number of chains in chain selection (" + to_string(ci.size()) +
            ") and base selection (" + to_string(cr.size()) + ")");
    }
    this->chainIndexes = ci;
    this->resIndexes = cr;
    this->withChainIndex = true;
    this->withResIndex = true;
}

BMSelection::BMSelection(const vector<int>& ci, const vector<vector<string> >& cr) {
    if(ci.size() != cr.size()) {
        throw invalid_argument("Inconssitent number of chains in chain selection (" + to_string(ci.size()) +
            ") and base selection (" + to_string(cr.size()) + ")");
    }
    this->chainIndexes = ci;
    this->chainResIDs = cr;
    this->withChainIndex = true;
    this->withResIndex = false;
}

BMSelection::BMSelection(const vector<string>& ci, const vector<vector<int> >& cr) {
    if(ci.size() != cr.size()) {
        throw invalid_argument("Inconssitent number of chains in chain selection (" + to_string(ci.size()) +
            ") and base selection (" + to_string(cr.size()) + ")");
    }
    this->chainIDs = ci;
    this->resIndexes = cr;
    this->withChainIndex = false;
    this->withResIndex = true;
}

BMSelection::BMSelection(const vector<string>& ci, const vector<vector<string> >& cr) {
    if(ci.size() != cr.size()) {
        throw invalid_argument("Inconssitent number of chains in chain selection (" + to_string(ci.size()) +
            ") and base selection (" + to_string(cr.size()) + ")");
    }
    this->chainIDs = ci;
    this->chainResIDs = cr;
    this->withChainIndex = false;
    this->withResIndex = false;
}

BMSelection::~BMSelection() {

}

BriqxModule::BriqxModule() {
    this->BMname = "NewModule";
}

BriqxModule::BriqxModule(const string& pdbFile, BasePairLib& bpl, AtomLib& atl,
    const string& name, const BMSelection& bms) {
    this->BMname = name;
    RNAPDB rnapdb(pdbFile, name);
    auto& chains_tmp = rnapdb.getChains();
    auto& bases_tmp = rnapdb.getBaseList();
    if(bms.getChainIndexes().size() == 0 && bms.getChainIDs().size() == 0) {
        // bms == NULL or empty, deep-copy whole RNAPDB, note the pdbID in RNAChains are not modified
        this->chains = move(chains_tmp);  // 使用vector移动构造，直接挪用 rnapdb 的资源
        this->baseList = move(bases_tmp);
        // string lastChainID = "@";
        // RNAChain* curChain;
        // int lb = bases_tmp.size();
        // map<string, RNAChain*> ID2chMap;
        // for(int i=0; i<lb;i++) {
        //     RNABase* iniBase = bases_tmp[i];
        //     RNABase* newBase = new RNABase(iniBase->baseID, iniBase->chainID, iniBase->baseType);
        //     newBase->setResSeqID(iniBase->baseSeqID);

        //     if(newBase->chainID != lastChainID) {
        //         try {
        //             curChain = ID2chMap.at(newBase->chainID);
        //         } catch(out_of_range) {
        //             curChain = new RNAChain(newBase->chainID);
        //             this->chains.emplace_back(curChain);
        //             ID2chMap.emplace(newBase->chainID, curChain);
        //         }
        //         lastChainID = newBase->chainID;
        //     }

        //     vector<Atom*>* iniAtomList = iniBase->getAtomList();
        //     int la = iniAtomList->size();
        //     for(int j=0;j<la;j++) {
        //         Atom* iniAtom = (*iniAtomList)[j];
        //         Atom* newAtom = new Atom();
        //         *newAtom = *iniAtom;
        //         newBase->addAtom(newAtom);
        //     }
        //     this->baseList.emplace_back(newBase);
        //     curChain->addBase(newBase);
        // }
    } else { // bms not empty, rebuild chains from bases
        bool isRsInd = bms.isUseResIndex();
        bool isChInd = bms.isUseChainIndex();  // 要对每一个找到的 RNABase 检查它的 chain 归属
        int lb = bases_tmp.size();
        map<int, RNABase*> ind2bsMap;
        map<int, int> ind2tmpMap;  // baseSeqID 到 base_tmp index 的 map
        map<string, RNABase*> ID2bsMap;
        map<string, int> ID2tmpMap;  // baseID 到 base_tmp index 的 map
        if(isRsInd) {
            for(int i=0; i<lb; i++) {
                ind2bsMap.emplace(bases_tmp[i]->baseSeqID, bases_tmp[i]);
                ind2tmpMap.emplace(bases_tmp[i]->baseSeqID, i);
            }
            const auto& bmsResInd = bms.getResIndexes();
            int l1 = bmsResInd.size();
            if(isChInd) { // chain index given
                const auto& bmsChInd = bms.getChainIndexes();
                for(int i=0;i<l1;i++) {
                    int l2 = bmsResInd[i].size();
                    string curChID;
                    try {
                        curChID = chains_tmp.at(bmsChInd.at(i))->getChainID();
                    } catch(out_of_range) {
                        throw out_of_range("chain index" + to_string(bmsChInd.at(i)) + " out of range");
                    }
                    RNAChain* curChain = new RNAChain(curChID);
                    for(int j=0;j<l2;j++) {
                        if(ind2bsMap.at(bmsResInd[i][j])->chainID == curChID) {
                            this->baseList.emplace_back(ind2bsMap.at(bmsResInd[i][j]));
                            bases_tmp[ind2tmpMap.at(bmsResInd[i][j])] = nullptr;
                            curChain->addBase(ind2bsMap.at(bmsResInd[i][j]));
                        } else {
                            throw invalid_argument("Base ("+ to_string(i) + "," + to_string(j) +
                                ") not in Chain indexed " + to_string(bmsChInd[i]));
                        }
                    }
                    this->chains.emplace_back(curChain);
                }
            } else { // chain ID given
                const auto& bmsChID = bms.getChainIDs();
                for(int i=0;i<l1;i++) {
                    int l2 = bmsResInd[i].size();
                    string curChID;
                    try {
                        curChID = bmsChID.at(i);
                    } catch(out_of_range) {
                        throw out_of_range("Speified chainID has less than "+to_string(i)+" chains");
                    }
                    RNAChain* curChain = new RNAChain(curChID);
                    for(int j=0;j<l2;j++) {
                        if(ind2bsMap.at(bmsResInd[i][j])->chainID == curChID) {
                            this->baseList.emplace_back(ind2bsMap.at(bmsResInd[i][j]));
                            bases_tmp[ind2tmpMap.at(bmsResInd[i][j])] = nullptr;
                            curChain->addBase(ind2bsMap.at(bmsResInd[i][j]));
                        } else {
                            throw invalid_argument("Base ("+ to_string(i) + "," + to_string(j) +
                                ") not in Chain ID " + curChID);
                        }
                    }
                    this->chains.emplace_back(curChain);
                }
            }
        } else { // isRsInd = false; residue ID given
            for(int i=0; i<lb; i++) {
                ID2bsMap.emplace(bases_tmp[i]->baseID, bases_tmp[i]);
                ID2tmpMap.emplace(bases_tmp[i]->baseID, i);
            }
            const auto& bmsResID = bms.getChainResIDs();
            int l1 = bmsResID.size();
            if(isChInd) { // chain index given
                const auto& bmsChInd = bms.getChainIndexes();
                for(int i=0;i<l1;i++) {
                    int l2 = bmsResID[i].size();
                    string curChID;
                    try {
                        curChID = chains_tmp.at(bmsChInd.at(i))->getChainID();
                    } catch(out_of_range) {
                        throw out_of_range("PDB has less than " + to_string(i) + " chains");
                    }
                    RNAChain* curChain = new RNAChain(curChID);
                    for(int j=0;j<l2;j++) {
                        if(ID2bsMap.at(bmsResID[i][j])->chainID == curChID) {
                            this->baseList.emplace_back(ID2bsMap.at(bmsResID[i][j]));
                            bases_tmp[ID2tmpMap.at(bmsResID[i][j])] = nullptr;
                            curChain->addBase(ID2bsMap.at(bmsResID[i][j]));
                        } else {
                            throw invalid_argument("Base ("+ to_string(i) + "," + to_string(j) +
                                ") not in Chain indexed " + to_string(bmsChInd[i]));
                        }
                    }
                    this->chains.emplace_back(curChain);
                }
            } else { // chain ID given
                const auto& bmsChID = bms.getChainIDs();
                for(int i=0;i<l1;i++) {
                    int l2 = bmsResID[i].size();
                    string curChID;
                    try {
                        curChID = bmsChID.at(i);
                    } catch(out_of_range) {
                        throw out_of_range("chains out of range");
                    }
                    RNAChain* curChain = new RNAChain(curChID);
                    for(int j=0;j<l2;j++) {
                        if(ID2bsMap.at(bmsResID[i][j])->chainID == curChID) {
                            this->baseList.emplace_back(ID2bsMap.at(bmsResID[i][j]));
                            bases_tmp[ID2tmpMap.at(bmsResID[i][j])] = nullptr;
                            curChain->addBase(ID2bsMap.at(bmsResID[i][j]));
                        } else {
                            throw invalid_argument("Base ("+ to_string(i) + "," + to_string(j) +
                                ") not in Chain ID " + curChID);
                        }
                    }
                    this->chains.emplace_back(curChain);
                }
            } // if isChInd
        } // if isResInd
    } // if bms empty
    int lb = this->baseList.size();
    for(int i=0;i<lb;i++) {
        for( int j=i+1;j<lb;j++) {
            if(baseList[i]->contactTo(baseList[j])) {
                this->basePairList.emplace_back(new BasePair(baseList[i], baseList[j], &atl));
            }
        }
    }
    int lbp = this->basePairList.size();
    for(int i=0;i<lbp;i++){
        BasePair* curBp = this->basePairList[i];
        BasePair* revBp = new BasePair(curBp->baseB, curBp->baseA, &atl);
        this->revBasePairList.emplace_back(revBp);
        this->bb2bpMap.emplace(array<RNABase*,2>{curBp->baseA, curBp->baseB}, curBp);
        this->bb2bpMap.emplace(array<RNABase*,2>{curBp->baseB, curBp->baseA}, revBp);
        this->PairEneMap.emplace(curBp, bpl.getPairEnergy(curBp->baseA, curBp->baseB));
        this->PairEneMap.emplace(revBp, bpl.getPairEnergy(revBp->baseA, revBp->baseB));
    }
}

BriqxModule::BriqxModule(const vector<RNAChain*>& chains, const vector<RNABase*>& baseList, BasePairLib& bpl,
    AtomLib& atl, const string& name, bool beLazy) {
    int lc = chains.size();
    for(int i=0;i<lc;i++) {
        this->chains.emplace_back(chains[i]);
    }
    int lb = baseList.size();
    for(int i=0;i<lb;i++) {
        this->baseList.emplace_back(baseList[i]);
    }

    if(beLazy) return;

    for(int i=0;i<lb;i++) {
        for( int j=i+1;j<lb;j++) {
            if(baseList[i]->contactTo(baseList[j])) {
                this->basePairList.emplace_back(new BasePair(baseList[i], baseList[j], &atl));
            }
        }
    }
    int lbp = this->basePairList.size();
    for(int i=0;i<lbp;i++){
        BasePair* curBp = this->basePairList[i];
        BasePair* revBp = new BasePair(curBp->baseB, curBp->baseA, &atl);
        this->revBasePairList.emplace_back(revBp);
        this->bb2bpMap.emplace(array<RNABase*,2>{curBp->baseA, curBp->baseB}, curBp);
        this->bb2bpMap.emplace(array<RNABase*,2>{curBp->baseB, curBp->baseA}, revBp);
        this->PairEneMap.emplace(curBp, bpl.getPairEnergy(curBp->baseA, curBp->baseB));
        this->PairEneMap.emplace(revBp, bpl.getPairEnergy(revBp->baseA, revBp->baseB));
    }
}

BriqxModule::BriqxModule(const BriqxModule& bm0, const TransForm& tf, const XYZ& Acog, const XYZ& Bcog,
    BasePairLib& bpl, AtomLib& atl, const string& name, bool beLazy)
{
    if(name == "") {
        BMname = bm0.BMname;
    } else {
        BMname = name;
    }
    string lastChainID = "@";
    RNAChain* curChain;
    int lb = bm0.baseList.size();
    map<string, RNAChain*> ID2chMap;
    for(int i=0; i<lb;i++) {
        RNABase* iniBase = bm0.baseList[i];
        RNABase* newBase = new RNABase(iniBase->baseID, iniBase->chainID, iniBase->baseType);
        newBase->setResSeqID(iniBase->baseSeqID);

        if(newBase->chainID != lastChainID) {
            try {
                curChain = ID2chMap.at(newBase->chainID);
            } catch(out_of_range) {
                curChain = new RNAChain(newBase->chainID);
                chains.emplace_back(curChain);
                ID2chMap.emplace(newBase->chainID, curChain);
            }
            lastChainID = newBase->chainID;
        }

        vector<Atom*>* iniAtomList = iniBase->getAtomList();
        int la = iniAtomList->size();
        for(int j=0;j<la;j++) {
            Atom* iniAtom = (*iniAtomList)[j];
            Atom* newAtom = new Atom();
            *newAtom = *iniAtom;
            newAtom->setCoord(tf.transform(newAtom->coord - Bcog) + Acog);
            newBase->addAtom(newAtom);
        }
        baseList.emplace_back(newBase);
        curChain->addBase(newBase);
    }

    if(beLazy) return;

    for(int i=0;i<lb;i++) {
        for( int j=i+1;j<lb;j++) {
            if(baseList[i]->contactTo(baseList[j])) {
                this->basePairList.emplace_back(new BasePair(baseList[i], baseList[j], &atl));
            }
        }
    }
    int lbp = this->basePairList.size();
    for(int i=0;i<lbp;i++){
        BasePair* curBp = this->basePairList[i];
        BasePair* revBp = new BasePair(curBp->baseB, curBp->baseA, &atl);
        this->revBasePairList.emplace_back(revBp);
        this->bb2bpMap.emplace(array<RNABase*,2>{curBp->baseA, curBp->baseB}, curBp);
        this->bb2bpMap.emplace(array<RNABase*,2>{curBp->baseB, curBp->baseA}, revBp);
        this->PairEneMap.emplace(curBp, bpl.getPairEnergy(curBp->baseA, curBp->baseB));
        this->PairEneMap.emplace(revBp, bpl.getPairEnergy(revBp->baseA, revBp->baseB));
    }
}

vector<BasePair*>& BriqxModule::getBasePairList(const BMSelection& bms, const bool withReversed) const {
    int lbp = this->basePairList.size();
    vector<BasePair*>* vec = new vector<BasePair*>;
    if(bms.isUseResIndex()) {
        map<int, set<BasePair*> > Ind2BpMap;
        for(int i=0;i<lbp;i++) {
            BasePair* curBp = this->basePairList[i];
            int IndA = curBp->baseA->baseSeqID;
            int IndB = curBp->baseB->baseSeqID;
            Ind2BpMap[IndA].emplace(curBp);
            Ind2BpMap[IndB].emplace(curBp);
        }
        if(withReversed) {
            for(int i=0;i<lbp;i++) {
                BasePair* curBp = this->revBasePairList[i];
                int IndA = curBp->baseA->baseSeqID;
                int IndB = curBp->baseB->baseSeqID;
                Ind2BpMap[IndA].emplace(curBp);
                Ind2BpMap[IndB].emplace(curBp);
            }
        }
        const vector<vector<int> >& bmsResInd = bms.getResIndexes();
        set<BasePair*> set0;
        int l1 = bmsResInd.size();
        for(int i=0; i<l1; i++) {
            int l2 = bmsResInd[i].size();
            for(int j=0;j<l2;j++) {
                set<BasePair*> curSet = Ind2BpMap.at(bmsResInd[i][j]);
                set0.insert(curSet.begin(), curSet.end()); // set0 is ordered by baseSeqID
            }
        }
        vec->assign(set0.begin(), set0.end());
    } else {
        map<string, set<BasePair*> > ID2BpMap;
        for(int i=0;i<lbp;i++) {
            BasePair* curBp = this->basePairList[i];
            string IDA = curBp->baseA->baseID;
            string IDB = curBp->baseB->baseID;
            ID2BpMap[IDA].emplace(curBp);
            ID2BpMap[IDB].emplace(curBp);
        }
        if(withReversed) {
            for(int i=0;i<lbp;i++) {
                BasePair* curBp = this->revBasePairList[i];
                string IDA = curBp->baseA->baseID;
                string IDB = curBp->baseB->baseID;
                ID2BpMap[IDA].emplace(curBp);
                ID2BpMap[IDB].emplace(curBp);
            }
        }
        const vector<vector<string> >& bmsResID = bms.getChainResIDs();
        set<BasePair*> set0;
        int l1 = bmsResID.size();
        for(int i=0; i<l1; i++) {
            int l2 = bmsResID[i].size();
            for(int j=0;j<l2;j++) {
                set<BasePair*> curSet = ID2BpMap.at(bmsResID[i][j]);
                set0.insert(curSet.begin(), curSet.end()); // set0 is ordered by baseSeqID
            }
        }
        vec->assign(set0.begin(), set0.end());
    }
    return *vec;
}

int BriqxModule::calcDDMMatrix(const BriqxModule& other, map<array<BasePair*, 2>, double>& DDMMatrix) const {
    const vector<BasePair*>& bpl1 = this->basePairList;
    const vector<BasePair*>& bplr1 = this->revBasePairList;
    const vector<BasePair*>& bpl2 = other.basePairList;
    int l1 = bpl1.size();
    int l2 = bpl2.size();
    for(int i=0;i<l1;i++) {
        for(int j=0;j<l2;j++) {
            BasePair* bp1 = bpl1[i];
            BasePair* bpr1 = bplr1[i];  // only reverse BasePair of THIS object
            BasePair* bp2 = bpl2[j];
            double DDM = bp1->dm.distanceTo(bp2->dm);
            array<BasePair*,2> key{bp1,bp2};
            DDMMatrix.emplace(key, DDM);
            double DDMr = bpr1->dm.distanceTo(bp2->dm);
            array<BasePair*,2> keyr{bpr1,bp2};
            DDMMatrix.emplace(keyr, DDMr);
        } // j loop
    } // i loop
    return EXIT_SUCCESS;
} // BriqxModule::calcDDMMatrix

int BriqxModule::sortMapByValue(const map<array<BasePair*, 2>, double>& unsortMap,
    vector<pair<array<BasePair*, 2>, double> >& sortedVec, string order) { //DONE: 改写成utils模板函数
    for(const auto& item : unsortMap) {
        sortedVec.emplace_back(item);
    }
    if(order == "desc") {
        sort(sortedVec.begin(), sortedVec.end(), [] (const auto& x, const auto& y) {return x.second > y.second;});
    } else if(order == "asc") {
        sort(sortedVec.begin(), sortedVec.end(), [] (const auto& x, const auto& y) {return x.second < y.second;});
    } else {
        // raise error
        throw invalid_argument("Unknown order " + order);
    }

    return EXIT_SUCCESS;
}

int BriqxModule::resolveTransformByAlign(const vector<array<BasePair*,2> >& alignVec, TransForm& tf, XYZ& Acog, XYZ& Bcog) {
        // first transfer BasePair to XYZ( 4-point representation for each base)
        vector<XYZ> points1;
        vector<XYZ> points2;
        int lv = alignVec.size();
	    if(lv == 0) return 0;
        for(int i=0;i<lv;i++) {
            const array<XYZ,4>& rep1A = alignVec[i][0]->baseA->getFourPseudoAtomCoords();  //常量左值引用可以用右值赋值
            points1.insert(points1.end(), rep1A.begin(), rep1A.end());
            const auto rep1B = alignVec[i][0]->baseB->getFourPseudoAtomCoords();
            points1.insert(points1.end(), rep1B.begin(), rep1B.end());
            const auto rep2A = alignVec[i][1]->baseA->getFourPseudoAtomCoords();
            points2.insert(points2.end(), rep2A.begin(), rep2A.end());
            const auto& rep2B = alignVec[i][1]->baseB->getFourPseudoAtomCoords();
            points2.insert(points2.end(), rep2B.begin(), rep2B.end());
        }
	    Acog = getCOG(points1);
	    Bcog = getCOG(points2);

	    vector<XYZ> listA;
	    vector<XYZ> listB;
        int lv8 = lv*8;
	    for(int i=0;i<lv8;i++){
	    	XYZ a = points1[i] - Acog;
	    	XYZ b = points2[i] - Bcog;
	    	listA.emplace_back(a);
	    	listB.emplace_back(b);
	    }

	    tf = buildRotation(listA, listB);
        return EXIT_SUCCESS;
}

int BriqxModule::resolveTransformByAlign(const vector<array<RNABase*,2> >& alignVec, TransForm& tf, XYZ& Acog, XYZ& Bcog) {
        // first transfer BasePair to XYZ( 4-point representation for each base)
        vector<XYZ> points1;
        vector<XYZ> points2;
        int lv = alignVec.size();
	    if(lv == 0) return 0;
        for(int i=0;i<lv;i++) {
            const auto& rep1A = alignVec[i][0]->getFourPseudoAtomCoords();
            points1.insert(points1.end(), rep1A.begin(), rep1A.end());
            const auto& rep2A = alignVec[i][1]->getFourPseudoAtomCoords();
            points2.insert(points2.end(), rep2A.begin(), rep2A.end());
        }
	    Acog = getCOG(points1);
	    Bcog = getCOG(points2);

	    vector<XYZ> listA;
	    vector<XYZ> listB;
        int lv4 = lv*4;
	    for(int i=0;i<lv4;i++){
	    	XYZ a = points1[i] - Acog;
	    	XYZ b = points2[i] - Bcog;
	    	listA.emplace_back(a);
	    	listB.emplace_back(b);
	    }

	    tf = buildRotation(listA, listB);
        return EXIT_SUCCESS;
}

int BriqxModule::coordTransform(const TransForm& tf, const XYZ& Acog,  const XYZ& Bcog, 
    vector<array<XYZ,4> >& tfBasePos) const
{
    tfBasePos.clear();
    int lb = this->baseList.size();
    if(lb==0) return 0;
    for(int i=0;i<lb;i++) {
        array<XYZ,4> iniCoord = this->baseList[i]->getFourPseudoAtomCoords();
        array<XYZ,4> ar0{};
        for(int j=0; j<4; j++) {
            ar0[j] = tf.transform(iniCoord[j]-Bcog) + Acog;
        }
        tfBasePos.emplace_back(ar0);
    }
    return EXIT_SUCCESS;
}

int BriqxModule::alignBasesByCoord(const BriqxModule& other, vector<array<RNABase*, 2> >& alignVec,
    const vector<array<XYZ,4> >& BasePos, const double& distMax) const {

    int lb1 = this->baseList.size();
    int lb2 = other.baseList.size();

    vector<array<XYZ, 4> > aPos;
    vector<array<XYZ, 4> > bPos;

    for(int i=0; i<lb1; i++) {
        aPos.emplace_back(baseList[i]->getFourPseudoAtomCoords());
    }
    int inBpSize = BasePos.size();
    if(inBpSize == 0) {
        // use native coordinates
        for(int i=0; i<lb2; i++) {
            bPos.emplace_back(other.baseList[i]->getFourPseudoAtomCoords());
        }
    } else if(inBpSize == lb2) {
        bPos = BasePos;
    } else {
        throw invalid_argument("invalid BasePos!");
    }

    vector<pair<pair<RNABase*,RNABase*>, double> > abmdVec{};
    map<RNABase*, RNABase*> abAlign{}, baAlign{};
    for(int i=0; i<lb1; i++) {
        for(int j=0;j<lb2;j++) {
            double dist = utils::meanDist<array<XYZ,4> >(aPos[i], bPos[j], 10, 1000);
            if(dist < distMax) {
                abmdVec.emplace_back(
                    pair<pair<RNABase*,RNABase*>, double>(
                        pair<RNABase*, RNABase*>(baseList[i], other.baseList[j]), dist
                    )
                );
            }
        }
    }
    sort(abmdVec.begin(), abmdVec.end(), [] (const auto& x, const auto& y) {return x.second < y.second;});
    int lv = abmdVec.size();
    for(int i=0;i<lv;i++) {
        pair<RNABase*, RNABase*> p1 = abmdVec[i].first;
        auto iter1 = abAlign.find(p1.first);
        auto iter2 = baAlign.find(p1.second);
        if(iter1 == abAlign.end() && iter2 == baAlign.end()) {
            abAlign.emplace(p1);
            baAlign.emplace(p1.second, p1.first);
        }
    }
    for(int i=0;i<lb1;i++) {
        RNABase* ba = baseList[i];
        auto iter1 = abAlign.find(ba);
        if(iter1 != abAlign.end()) {
            alignVec.emplace_back(array<RNABase*,2>{ba, iter1->second});
        }
    }
    return EXIT_SUCCESS;
}

void BriqxModule::printPDBFormat(ofstream& out) const {
    int startID = 1;
    int lc = chains.size();
    for(unsigned int i=0;i<lc;i++)
    {
    	RNAChain* pc = this->chains.at(i);
        int lb = pc->getChainLength();
        auto& bl = pc->getBaseList();
        for(int j=0;j<lb;j++)
        {
            RNABase* res = bl.at(j);
            startID = res->printPDBFormat(out,startID);
        }
    }
}

void BriqxModule::deepClear() {
    int lb = baseList.size();
    for(int i=0;i<lb;i++) {
        delete baseList[i];
    }
    int lbp = basePairList.size();
    for(int i=0;i<lbp;i++) {
        delete basePairList[i];
    }
    int lrbp = revBasePairList.size();
    for(int i=0;i<lrbp;i++) {
        delete revBasePairList[i];
    }
    int lc = chains.size();
    for(int i=0;i<lc;i++) {
        delete chains[i];
    }
}

void BriqxModule::clearBasePairs() {
    int lbp = basePairList.size();
    for(int i=0;i<lbp;i++) {
        delete basePairList[i];
    }
    int lrbp = revBasePairList.size();
    for(int i=0;i<lrbp;i++) {
        delete revBasePairList[i];
    }
}

BriqxModule::~BriqxModule() {
    //TODO: do nothing
}
}  // namespace NSPbm
