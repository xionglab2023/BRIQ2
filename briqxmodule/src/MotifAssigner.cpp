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
            motifInfo = new graphInfo(nNode, seq, connectToDownstream, allNodes, ene, fullGraph->atLib);
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


    void MotifAssigner::bySeed() {
        vector<NuEdge*> seedEdge;  // take non-neighboring nonWC pairs as seed
        int seqLen = nuGraph->seqLen;
        for(int i=0; i<seqLen; i++) {
            for(int j=i+2; j<seqLen; j++) {
                NuEdge* curEdge = nuGraph->allEdges[i*seqLen+j];
                NuEdge* revEdge = nuGraph->allEdges[j*seqLen+i];

                // Exclude WC Edge in seed search
                if(curEdge->isWC()) {
                    continue;
                }
                
                if(curEdge->weight <= revEdge->weight && curEdge->weight <= MOTIFASSIGNER_SEED_CUTOFF) {
                    seedEdge.emplace_back(curEdge);
                } else if(revEdge->weight < curEdge->weight && revEdge->weight <= MOTIFASSIGNER_SEED_CUTOFF) {
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
        for(int i=0; i<nSeed; i++) {
            set<int> selNode{seedEdge[i]->indexA, seedEdge[i]->indexB};
            set<int>* pCurConnect = new set<int>;
            set<int>* pLastConnect = new set<int>, *ptmp;
            set<NuEdge*> selEdge = {seedEdge[i]};
            for(auto& it:selNode) {
                // for(const int& it2: *(nnbNeighborMap[it])) {
                //     connectNode.emplace(it2);
                // }
                set_union(pCurConnect->begin(), pCurConnect->end(), 
                          nnbNeighborMap[it]->begin(), nnbNeighborMap[it]->end(),
                          inserter(*pCurConnect, pCurConnect->begin()));
                for(const int& it2: *(nbNeighborMap[it])) {
                    if(WCNodeSet.count(it2) && WCNodeSet.count(it)) {
                        continue;
                    }
                    pCurConnect->emplace(it2);
                }
            }
            while(!includes(selNode.begin(), selNode.end(), pCurConnect->begin(), pCurConnect->end())) {
                set_union(selNode.begin(), selNode.end(),
                          pCurConnect->begin(), pCurConnect->end(),
                          inserter(selNode, selNode.begin()));
                //TODO selEdge
                pLastConnect->clear();
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
                ptmp = pCurConnect;
                pCurConnect = pLastConnect;
                pLastConnect = ptmp;
            }
            delete pCurConnect;
            delete pLastConnect;
            
            // TODO: filter motifs by recheck compactness

            vector<int> NodeIndexVec(selNode.begin(), selNode.end());

            // TODO: compute ene
            double ene = 0;

            MotifGraph* pMotif = new MotifGraph(nuGraph, NodeIndexVec, ene);

            // TODO: redundent reduction
            int nmtf =  motifs.size();
            bool dup = false;
            for(int i=0; i<nmtf; i++) {
                if(*motifs[i] == *pMotif) {
                    dup = true;
                    break;
                }
            }
            dup || motifs.emplace_back(pMotif);
            
        }

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