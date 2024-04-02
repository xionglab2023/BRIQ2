/**
 * @file MotifAssigner.cpp
 * @author Kalrk Chen (klarkchen@ustc.edu.cn)
 * @brief  Implementation of MotifAssigner class
 * @version 0.1
 * @date 2024-01-03
 * 
 * @copyright Copyright (c) 2024 XLAB
 * 
 */

#include "briqxmodule/MotifAssigner.h"
#include <string.h>
#include <algorithm>

namespace NSPbm {
    using namespace NSPpredNA;
    using namespace std;

    
    MotifGraph::MotifGraph(NuGraph* graphIn, vector<int>& nodeIndexVec, double eneIn, set<NuEdge*>* pHelixEdge) {
        fullGraph = graphIn;
        nNode = nodeIndexVec.size();
        nFrag = 1;
        int iHelix = 0;
        sep.emplace_back(0);
        ene = eneIn;
        seq = new int[nNode];
        nodeIndex = new int[nNode];
        memcpy(nodeIndex, &nodeIndexVec[0], sizeof(int)*nNode);
        connectToDownstream = new bool[nNode];
        fixed = new bool[nNode];
        allNodes = new NuNode*[nNode];
        allEdges = new NuEdge*[nNode*nNode];
        geList.reserve(nNode*(nNode-1)/2);
        set<int> NodeIndexSet(nodeIndexVec.begin(), nodeIndexVec.end());
        map<int,NuEdge*>* pHelixNode2EdgeMap = new map<int,NuEdge*>;
        for(auto& it: *pHelixEdge) {
            if(pHelixNode2EdgeMap->count(it->indexA)) {
                string errInfo = "[ERROR] Duplicated Helix Node " + to_string(it->indexA) + " in edge" +\
                    to_string(it->indexA) + "-" + to_string(it->indexB) + "\n";
                cerr << errInfo;
            }
            if(pHelixNode2EdgeMap->count(it->indexB)) {
                string errInfo = "[ERROR] Duplicated Helix Node " + to_string(it->indexB) + " in edge" +\
                    to_string(it->indexA) + "-" + to_string(it->indexB) + "\n";
                cerr << errInfo;
            }
            pHelixNode2EdgeMap->emplace(it->indexA, it);
            pHelixNode2EdgeMap->emplace(it->indexB, it);
        }
        map<int,int> node2HelixMap;  // map from nodeIndex to ihelix
        for(int i=0; i<nNode; i++) {
            allNodes[i] = graphIn->allNodes[nodeIndex[i]];
            seq[i] = graphIn->seq[nodeIndex[i]];
            if(NodeIndexSet.count(nodeIndex[i]+1)) {
                connectToDownstream[i] = graphIn->connectToDownstream[nodeIndex[i]];
                if(!connectToDownstream[i]) {
                    sep.emplace_back(0);
                    sep.emplace_back(0);
                    nFrag+=1;
                }
            } else if(nodeIndex[i]+1 == graphIn->seqLen) {
                connectToDownstream[i] = false;
                sep.emplace_back(0);
            }
            else {
                connectToDownstream[i] = false;  // the original downstream node is not included in subgraph
                // next motif node is nodeIndexVec[i+1], next full seq node is nodeIndexVec[i]+1
                // check helix first
                if(pHelixNode2EdgeMap->count(nodeIndexVec[i]) && pHelixNode2EdgeMap->count(nodeIndexVec[i]+1)
                    && graphIn->allEdges[nodeIndexVec[i]*(graphIn->seqLen)+nodeIndexVec[i]+1]->weight
                    < MOTIFASSIGNER_GROW_CUTOFF) {  // helix break
                    iHelix+=1;
                    sep.emplace_back(-iHelix);
                    node2HelixMap.emplace(nodeIndexVec[i], iHelix);
                }


            }
            fixed[i] = graphIn->fixed[nodeIndex[i]];
            for(int j=0; j<nNode; j++) {
                allEdges[i*nNode+j] = graphIn->allEdges[nodeIndex[i]*graphIn->seqLen+nodeIndex[j]];
            }
            for(int j=i+1; j<nNode; j++) {
                geList.emplace_back(allEdges[i*nNode+j]);
            }
        }
        motifInfo = nullptr;
    }

    int MotifGraph::writePDB(const string& outputFile) {
        if(!motifInfo) {
            motifInfo = new graphInfo(nNode, seq, connectToDownstream, fixed, allNodes, ene, fullGraph->atLib, 0);
        }
        motifInfo->printPDB(outputFile);
        return EXIT_SUCCESS;
    }


    bool MotifGraph::operator==(const MotifGraph& other) {
        if(nNode != other.nNode || fullGraph != other.fullGraph) {
            return false;
        }
        for(int i=0; i<nNode; i++) {
            if(nodeIndex[i] != other.nodeIndex[i]) {
                return false;
            }
        }
        return true;
    }


    MotifAssigner::MotifAssigner(NuGraph* nuGraphIn) {
        nuGraph = nuGraphIn;
        int seqLen = nuGraph->seqLen;
        baseStackingMap = getStackingSetMap();
        for(int i=0; i<seqLen; i++) {
            for(int j=i+2; j<seqLen; j++) {
                NuEdge* curEdge = nuGraph->allEdges[i*seqLen+j];
                NuEdge* revEdge = nuGraph->allEdges[j*seqLen+i];
                if(isHelixByWC(i,j)) {
                    helixEdge.emplace(curEdge);
                    helixEdge.emplace(revEdge);
                }
                // helix by stacking
                /*
                if(isHelixByStacking(i,j)) {
                    helixEdge.emplace(curEdge);
                    helixEdge.emplace(revEdge);
                } 
                */
            }
        }
    }

    bool MotifAssigner::isHelixByWC(int i, int j) {
        int seqLen = nuGraph->seqLen;
        bool isWCIJ = false;
        bool isWCILast1JNext1 =false, isWCINext1JLast1=false;
        bool isWCILast2JNext2 =false, isWCINext2JLast2=false;
        if(i>j) {
            int tmp = j;
            i=j;
            j=tmp;
        }
        int seqLen = nuGraph->seqLen;
        isWCIJ = nuGraph->allEdges[i*seqLen + j]->isWC();
        if(isWCIJ) {
            if(i>1 && j<seqLen-2 && nuGraph->connectToDownstream[i-1] && nuGraph->connectToDownstream[i-2] &&
            nuGraph->connectToDownstream[j] && nuGraph->connectToDownstream[j+1]) {
                isWCILast2JNext2 =  nuGraph->allEdges[(i-2)*seqLen+j+2]->isWC();
                isWCILast1JNext1 =  nuGraph->allEdges[(i-1)*seqLen+j+1]->isWC();
            } else if (i==1 && j<seqLen-1 && nuGraph->connectToDownstream[i-1] && nuGraph->connectToDownstream[j]) {
                isWCILast1JNext1 = nuGraph->allEdges[(i-1)*seqLen+j+1]->isWC();
            }

            if(i<seqLen-2 && nuGraph->connectToDownstream[i] && nuGraph->connectToDownstream[i+1] && 
            nuGraph->connectToDownstream[j-1] && nuGraph->connectToDownstream[j-2]) {
                isWCINext2JLast2 = nuGraph->allEdges[(i+2)*seqLen+j-2]->isWC();
                isWCINext1JLast1 = nuGraph->allEdges[(i+1)*seqLen+j-1]->isWC();
            } else if(i == seqLen-2 && nuGraph->connectToDownstream[i] && nuGraph->connectToDownstream[j-1]) {
                isWCINext1JLast1 = nuGraph->allEdges[(i+1)*seqLen+j-1]->isWC();
            }
            if(isWCINext1JLast1 && isWCILast1JNext1) {
                return true;
            }
            if(isWCINext2JLast2 && isWCINext1JLast1) {
                return true;
            }
            if(isWCILast2JNext2 && isWCILast1JNext1) {
                return true;
            }
        }
        return false;
    }


    map<int, set<int>*>* MotifAssigner::getStackingSetMap() {
        int seqLen = nuGraph->seqLen;
        AtomLib* atl = nuGraph->atLib;
        vector<RNABase*> baseVec;
        map<int,char> baseTypeitoaMap = {{0,'A'}, {1, 'U'}, {2, 'G'}, {3, 'C'},
            {4,'a'}, {5, 't'}, {6, 'g'}, {7,'c'}, {-1, 'N'}};
        for(int i=0; i<seqLen;i++) {
            char btc = baseTypeitoaMap[nuGraph->allNodes[i]->baseType];
            RNABase* curBase = new RNABase(to_string(i),"A",btc);
            vector<Atom*> atomList = nuGraph->allNodes[i]->toAtomListWithPho(atl);
            int nAt = atomList.size();
            for(int j=0; j<nAt; j++) {
                curBase->addAtom(atomList[j]);
            }
            baseVec.emplace_back(curBase);
        }
        map<int, set<int>*> * baseStackingMap = new map<int, set<int>*>;
        for(int i=0; i<seqLen; i++) {
            for(int j=i+1; j<seqLen; j++) {
                if(baseVec[i]->isStackingTo(baseVec[j], atl)) {
                    if(baseStackingMap->count(i)==0) {
                        baseStackingMap->emplace(i,new set{j});
                    } else {
                        baseStackingMap->at(i)->emplace(j);
                    }
                    if(baseStackingMap->count(j)==0) {
                        baseStackingMap->emplace(j,new set{i});
                    } else {
                        baseStackingMap->at(j)->emplace(i);
                    }
                }
            }
        }

        // release baseVec pointed memory
        for(int i=0; i<seqLen;i++) {
            RNABase* p1 = baseVec[i];
            int nAt = p1->getAtomList()->size();
            for(int j=0; j<nAt; j++) {
                delete p1->getAtomList()->at(j);
            }
            delete p1;
        }

        return baseStackingMap;

    }


    bool MotifAssigner::isHelixByStacking(int i, int j) {
        int seqLen = nuGraph->seqLen;
        bool isStackILast1=false, isStackINext1=false, isStackJLast1=false, isStackJNext1=false;
        bool isStackILast2=false, isStackINext2=false, isStackJLast2=false, isStackJNext2=false;
        bool isBpILast1JNext1 =false, isBpINext1JLast1=false;
        bool isBpILast2JNext2 =false, isBpINext2JLast2=false;
        if(nuGraph->allEdges[i*seqLen+j]->weight < MOTIFASSIGNER_BP_WEIGHT_MAX) {
            if(i>0 && nuGraph->connectToDownstream[i-1]) {
                if(baseStackingMap->count(i) && baseStackingMap->at(i)->count(i-1)) {
                    isStackILast1 = true;
                    if(i>1 && nuGraph->connectToDownstream[i-2]) {
                        if(baseStackingMap->at(i-1)->count(i-2)) {
                            isStackILast2 = true;
                        }
                    }
                }
            }
            if(j>0 && nuGraph->connectToDownstream[j-1]) {
                if(baseStackingMap->count(j) && baseStackingMap->at(j)->count(j-1)) {
                    isStackJLast1 = true;
                    if(j>1 && nuGraph->connectToDownstream[j-2]) {
                        if(baseStackingMap->at(j-1)->count(j-2)) {
                            isStackJLast2 = true;
                        }
                    }
                }
            }
            if(nuGraph->connectToDownstream[i]) {
                if(baseStackingMap->count(i+1) && baseStackingMap->at(i+1)->count(i)) {
                    isStackINext1 = true;
                    if(i<seqLen-2 && nuGraph->connectToDownstream[i+1]) {
                        if(baseStackingMap->at(i+2)->count(i+1)) {
                            isStackINext2 = true;
                        }
                    }
                }
            }
            if(nuGraph->connectToDownstream[j]) {
                if(baseStackingMap->count(j+1) && baseStackingMap->at(j+1)->count(j)) {
                    isStackJNext1 = true;
                    if(j<seqLen-2 && nuGraph->connectToDownstream[j+1]) {
                        if(baseStackingMap->at(j+2)->count(j+1)) {
                            isStackJNext2 = true;
                        }
                    }
                }
            }

            if(isStackILast1 && isStackJNext1) {
                if(nuGraph->allEdges[(i-1)*seqLen+j+1]->weight < MOTIFASSIGNER_BP_WEIGHT_MAX &&
                   nuGraph->allEdges[(j+1)*seqLen+i-1]->weight < MOTIFASSIGNER_BP_WEIGHT_MAX) {
                    isBpILast1JNext1=true;
                }
            }
            if(isStackINext1 && isStackJLast1) {
                if(nuGraph->allEdges[(i+1)*seqLen+j-1]->weight < MOTIFASSIGNER_BP_WEIGHT_MAX &&
                   nuGraph->allEdges[(j-1)*seqLen+i+1]->weight < MOTIFASSIGNER_BP_WEIGHT_MAX) {
                    isBpINext1JLast1=true;
                }
            }
        
            if(isBpILast1JNext1 && isBpINext1JLast1) {
                return true;
            }

            if(isStackILast2 && isStackJNext2) {
                if(nuGraph->allEdges[(i-2)*seqLen+j+2]->weight < MOTIFASSIGNER_BP_WEIGHT_MAX &&
                   nuGraph->allEdges[(j+2)*seqLen+i-2]->weight < MOTIFASSIGNER_BP_WEIGHT_MAX) {
                    isBpILast2JNext2=true;
                }
            }

            if(isBpILast1JNext1 && isBpILast2JNext2) {
                return true;
            }

            if(isStackINext2 && isStackJLast2) {
                if(nuGraph->allEdges[(i+2)*seqLen+j-2]->weight < MOTIFASSIGNER_BP_WEIGHT_MAX &&
                   nuGraph->allEdges[(j-2)*seqLen+i+2]->weight < MOTIFASSIGNER_BP_WEIGHT_MAX) {
                    isBpINext2JLast2=true;
                }
            }

            if(isBpINext1JLast1 && isBpINext2JLast2) {
                return true;
            }
        }
        return false;
    }

    /**
     * @brief Determine if edge between node i,j is formed due to implicit interaction. An edge is
     * thought implicit if node i/j forms basepair with a third node k that stacking to j/i with
     * a stronger basepair weight than edge ij.
     * 
     * @param i 
     * @param j 
     * @return true 
     * @return false 
     */
    bool MotifAssigner::isImplicitEdge(int i, int j) {
        int seqLen = nuGraph->seqLen;
        double curWeight = min(nuGraph->allEdges[i*seqLen+j]->weight,
                               nuGraph->allEdges[j*seqLen+i]->weight);
        if(baseStackingMap->count(i)) {
            for(auto& it : *(baseStackingMap->at(i))) {
                if(min(nuGraph->allEdges[it*seqLen+j]->weight,
                       nuGraph->allEdges[j*seqLen+it]->weight) < curWeight) {
                    return true;
                }
            }
        }
        if(baseStackingMap->count(j)) {
            for(auto& it : *(baseStackingMap->at(j))) {
                if(min(nuGraph->allEdges[it*seqLen+i]->weight,
                       nuGraph->allEdges[i*seqLen+it]->weight) < curWeight) {
                    return true;
                }
            }
        }
        return false;
    }

    int MotifAssigner::writeEdgeWeight(ostream& outCSV) {
        int seqLen = nuGraph->seqLen;
        outCSV << "nodeID1,nodeID2,Weight,isWC,isNeighbor" << endl;
        for(int i=0; i<seqLen; i++) {
            for(int j=i+1; j<seqLen; j++) {
                NuEdge* curEdge = nuGraph->allEdges[i*seqLen+j];
                NuEdge* revEdge = nuGraph->allEdges[j*seqLen+i];
                outCSV << i << "," 
                       << j << ","
                       << curEdge->weight << ","
                       << noboolalpha << curEdge->isWC() << ","
                       << noboolalpha << (abs(i-j) == 1) << '\n';
                outCSV << j << "," 
                       << i << ","
                       << revEdge->weight << ","
                       << noboolalpha << revEdge->isWC() << ","
                       << noboolalpha << (abs(i-j) == 1) << '\n';
            }
        }
        return EXIT_SUCCESS;
    }

    #ifdef DEBUG
    void MotifAssigner::bySeed(string&& outPath, string&& outPDBprefix) {
    #else
    void MotifAssigner::bySeed() {
    #endif
        vector<NuEdge*> seedEdge;  // take non-neighboring nonWC pairs as seed
        int seqLen = nuGraph->seqLen;
        for(int i=0; i<seqLen; i++) {
            for(int j=i+2; j<seqLen; j++) {
                NuEdge* curEdge = nuGraph->allEdges[i*seqLen+j];
                NuEdge* revEdge = nuGraph->allEdges[j*seqLen+i];

                // Exclude helix WC Edge in seed search
                if(curEdge->isWC() && helixEdge.count(curEdge)) {
                    continue;
                }
                
                if(curEdge->weight <= revEdge->weight && curEdge->weight <= MOTIFASSIGNER_SEED_CUTOFF) {
                    cout << "get seed " << i << "-" << j << ", weight=" << curEdge->weight << ", WC: " << curEdge->isWC() << endl;
                    seedEdge.emplace_back(curEdge);
                } else if(revEdge->weight < curEdge->weight && revEdge->weight <= MOTIFASSIGNER_SEED_CUTOFF) {
                    cout << "get seed " << j << "-" << i << ", weight=" << revEdge->weight << ", WC: " << revEdge->isWC() << endl;
                    seedEdge.emplace_back(revEdge);
                }
            }
        }
        sort(seedEdge.begin(), seedEdge.end(),
            [](const auto& x, const auto& y){return x->weight < y->weight;});
        
        set<NuNode*> allNodeSet(nuGraph->allNodes, nuGraph->allNodes+seqLen);
        set<NuEdge*> allEdgeSet(nuGraph->allEdges, nuGraph->allEdges+seqLen*seqLen);
        map<int, set<NuEdge*>*> outEdgeMap, inEdgeMap;  // map from Node index to Edge vector
        for(int i=0;i<seqLen;i++) {
            if(!outEdgeMap.count(i)) {
                outEdgeMap.emplace(i, new set<NuEdge*>);
            }
            for(int j=0;j<seqLen;j++) {
                if(i==j) continue;
                if(!inEdgeMap.count(j)) {
                    inEdgeMap.emplace(j, new set<NuEdge*>);
                }
                if(nuGraph->allEdges[i*seqLen+j]->weight < MOTIFASSIGNER_GROW_CUTOFF) {
                    if(isImplicitEdge(i,j)) continue;
                    outEdgeMap[i]->emplace(nuGraph->allEdges[i*seqLen+j]);
                    inEdgeMap[j]->emplace(nuGraph->allEdges[i*seqLen+j]);
                }
            }
        }

        // map from nodes to the set of nodes that connect to the key node by an edge
        map<int, set<int>*> nbNeighborMap, nnbNeighborMap; 
        for(int i=0;i<seqLen;i++) {
            nbNeighborMap.emplace(i,new set<int>);
            nnbNeighborMap.emplace(i,new set<int>);
            for(auto& it : *outEdgeMap[i]) {
                if(it->sep < 2) {
                    nbNeighborMap[i]->emplace(it->indexB);
                } else {
                    nnbNeighborMap[i]->emplace(it->indexB);
                }
            }
            for(auto& it : *inEdgeMap[i]) {
                if(it->sep < 2) {
                    nbNeighborMap[i]->emplace(it->indexA);
                } else {
                    nnbNeighborMap[i]->emplace(it->indexA);
                }
            }
        }

        set<int> WCNodeSet;  // Set of nodes that have only WC non-neighbor basepair
        for(auto& it: nnbNeighborMap) {
            if(it.second->size() == 1 && 
            nuGraph->allEdges[it.first*seqLen + *(it.second->begin())]->isWC()) {
                WCNodeSet.emplace(it.first);
            }
        }

        int nSeed = seedEdge.size();
        #ifdef DEBUG
            ofstream outCSV;
            string csvFileName = outPath + "/" + outPDBprefix + "-Seeds" + ".csv";
            outCSV.open(csvFileName, ios::out);
            if(!outCSV.is_open()) {
                throw "[Error] Fail to open file" + csvFileName;
            }
            outCSV << "indexA,indexB,weight,isWC" << endl;
        #endif
        for(int i=0; i<nSeed; i++) {
            set<int> selNode{seedEdge[i]->indexA, seedEdge[i]->indexB};
            #ifdef DEBUG
                cout << "seed node: " << seedEdge[i]->indexA << "," << seedEdge[i]->indexB << endl;
                outCSV << seedEdge[i]->indexA << ","
                       << seedEdge[i]->indexB << ","
                       << seedEdge[i]->weight << ","
                       << seedEdge[i]->isWC() << endl;
            #endif
            set<int>* pCurConnect = new set<int>;
            set<int>* pLastConnect = new set<int>, *ptmp;
            set<NuEdge*> selEdge = {seedEdge[i]};
            // Initialize pCurConnect set
            for(auto& it:selNode) {
                // for(const int& it2: *(nnbNeighborMap[it])) {
                //     connectNode.emplace(it2);
                // }
                // connect all nnb neighbors (non-sequential-neighboring interactions)
                set_union(pCurConnect->begin(), pCurConnect->end(), 
                          nnbNeighborMap[it]->begin(), nnbNeighborMap[it]->end(),
                          inserter(*pCurConnect, pCurConnect->begin()));
                // connect nb neighbors
                for(const int& it2: *(nbNeighborMap[it])) {
                    // if current node is WC-only node and the nb neighbor is also WC-only node, then grow will stop,
                    // so the nb neighbor excluded
                    if(WCNodeSet.count(it2) && WCNodeSet.count(it)) {
                        continue;
                    }
                    pCurConnect->emplace(it2);
                }
            }
            // iterate until pCurConnect nodes are already included in selNode
            while(!includes(selNode.begin(), selNode.end(), pCurConnect->begin(), pCurConnect->end())) {
                set_union(selNode.begin(), selNode.end(),
                          pCurConnect->begin(), pCurConnect->end(),
                          inserter(selNode, selNode.begin()));
                //TODO selEdge
                pLastConnect->clear();
                // initialize the connected nodes in next round according to pCurConnect, use pLastConnect temporarily
                for(auto& it:*pCurConnect) {
                    // for(const int& it2: *(nnbNeighborMap[it])) {
                    //     connectNode.emplace(it2);
                    // }
                    set_union(pLastConnect->begin(), pLastConnect->end(), 
                              nnbNeighborMap[it]->begin(), nnbNeighborMap[it]->end(),
                              inserter(*pLastConnect, pLastConnect->begin()));
                    for(const int& it2: *(nbNeighborMap[it])) {
                        if(WCNodeSet.count(it2) && WCNodeSet.count(it)) {
                            continue;
                        }
                        pLastConnect->emplace(it2);
                    }
                }
                // point pLastConnect to the last round pCurConnect, and pCurConnect to the next round connected nodes
                ptmp = pCurConnect;
                pCurConnect = pLastConnect;
                pLastConnect = ptmp;
            }
            delete pCurConnect;
            delete pLastConnect;
            
            // TODO: filter motifs by recheck compactness

            vector<int> nodeIndexVec(selNode.begin(), selNode.end()); // set is internally sorted
            #ifdef DEBUG
                int nSelNode = nodeIndexVec.size();
                cout << "Selected Nodes: ";
                for(int j=0; j<nSelNode; j++) {
                    cout << nodeIndexVec[j] << ",";
                }
                cout << endl;
            #endif

            // TODO: compute ene
            double ene = 0;

            MotifGraph* pMotif = new MotifGraph(nuGraph, nodeIndexVec, ene);

            // TODO: redundent reduction
            int nmtf =  motifs.size();
            bool dup = false;
            for(int j=0; j<nmtf; j++) {
                if(*motifs[j] == *pMotif) {
                    dup = true;
                    break;
                }
            }
            dup || motifs.emplace_back(pMotif);
            
        }
        #ifdef DEBUG
            outCSV.close();
        #endif

        for(auto& it:outEdgeMap) {
            delete it.second;
        }
        for(auto& it:inEdgeMap) {
            delete it.second;
        }
        for(auto& it:nbNeighborMap) {
            delete it.second;
        }
        for(auto& it:nnbNeighborMap) {
            delete it.second;
        }
    }

    void MotifAssigner::byLouvain() {
        return;
    }

    MotifAssigner::~MotifAssigner() {
        //TODO
        for(auto& it : motifs) {
            delete it;
        }
        for(auto& it : *baseStackingMap) {
            delete it.second;
        }
        delete baseStackingMap;
    }
}