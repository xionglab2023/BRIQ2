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

    MotifGraph::MotifGraph(NuGraph* graphIn, vector<int>& NodeIndexVec, double eneIn) {
        fullGraph = graphIn;
        nNode = NodeIndexVec.size();
        ene = eneIn;
        seq = new int[nNode];
        nodeIndex = new int[nNode];
        memcpy(nodeIndex, &NodeIndexVec[0], sizeof(int)*nNode);
        connectToDownstream = new bool[nNode];
        fixed = new bool[nNode];
        allNodes = new NuNode*[nNode];
        allEdges = new NuEdge*[nNode*nNode];
        geList.reserve(nNode*(nNode-1)/2);
        set<int> NodeIndexSet(NodeIndexVec.begin(), NodeIndexVec.end());
        for(int i=0; i<nNode; i++) {
            allNodes[i] = graphIn->allNodes[nodeIndex[i]];
            seq[i] = graphIn->seq[nodeIndex[i]];
            if(NodeIndexSet.count(nodeIndex[i]+1)) {
                connectToDownstream[i] = graphIn->connectToDownstream[nodeIndex[i]];
            } else {
                connectToDownstream[i] = false;  // the original downstream node is not included in subgraph
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
        for(int i=0; i<seqLen; i++) {
            for(int j=i+2; j<seqLen; j++) {
                NuEdge* curEdge = nuGraph->allEdges[i*seqLen+j];
                NuEdge* revEdge = nuGraph->allEdges[j*seqLen+i];
                bool isStackILast=false, isStackINext=false, isStackJLast=false, isStackJNext=false;
                bool isBpILastJNext =false, isBpINextJLast=false;
                if(i>0 && nuGraph->connectToDownstream[i-1]) {
                    if(nuGraph->allEdges[(i-1)*seqLen+i]->weight < MOTIFASSIGNER_STACK_WEIGHT_MAX &&
                       nuGraph->allEdges[i*seqLen+i-1]->weight < MOTIFASSIGNER_STACK_WEIGHT_MAX) {
                        isStackILast = true;
                    }
                }
                if(nuGraph->connectToDownstream[j-1]) {
                    if(nuGraph->allEdges[(j-1)*seqLen+j]->weight < MOTIFASSIGNER_STACK_WEIGHT_MAX &&
                       nuGraph->allEdges[j*seqLen+j-1]->weight < MOTIFASSIGNER_STACK_WEIGHT_MAX) {
                        isStackJLast = true;
                    }
                }
                if(nuGraph->connectToDownstream[i]) {
                    if(nuGraph->allEdges[(i+1)*seqLen+i]->weight < MOTIFASSIGNER_STACK_WEIGHT_MAX &&
                       nuGraph->allEdges[i*seqLen+i+1]->weight < MOTIFASSIGNER_STACK_WEIGHT_MAX) {
                        isStackINext = true;
                    }

                }
                if(nuGraph->connectToDownstream[j]) {
                    if(nuGraph->allEdges[(j+1)*seqLen+j]->weight < MOTIFASSIGNER_STACK_WEIGHT_MAX &&
                       nuGraph->allEdges[j*seqLen+j+1]->weight < MOTIFASSIGNER_STACK_WEIGHT_MAX) {
                        isStackJNext = true;
                    }
                }
                if(isStackILast && isStackJNext) {
                    if(nuGraph->allEdges[(i-1)*seqLen+j+1]->weight < MOTIFASSIGNER_BP_WEIGHT_MAX &&
                       nuGraph->allEdges[(j+1)*seqLen+i-1]->weight < MOTIFASSIGNER_BP_WEIGHT_MAX) {
                        isBpILastJNext=true;
                    }
                }
                if(isStackINext && isStackJLast) {
                    if(nuGraph->allEdges[(i+1)*seqLen+j-1]->weight < MOTIFASSIGNER_BP_WEIGHT_MAX &&
                       nuGraph->allEdges[(j-1)*seqLen+i+1]->weight < MOTIFASSIGNER_BP_WEIGHT_MAX) {
                        isBpINextJLast=true;
                    }
                }
                
                if((isStackILast && isStackJNext && isBpILastJNext) || 
                   (isStackINext && isStackJLast && isBpINextJLast)) {
                    helixEdge.emplace(curEdge);
                    helixEdge.emplace(revEdge);
                }
            }
        }
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

            vector<int> NodeIndexVec(selNode.begin(), selNode.end());
            #ifdef DEBUG
                int nSelNode = NodeIndexVec.size();
                cout << "Selected Nodes: ";
                for(int j=0; j<nSelNode; j++) {
                    cout << NodeIndexVec[j] << ",";
                }
                cout << endl;
            #endif

            // TODO: compute ene
            double ene = 0;

            MotifGraph* pMotif = new MotifGraph(nuGraph, NodeIndexVec, ene);

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
    }
}