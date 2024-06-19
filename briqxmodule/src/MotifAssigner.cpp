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
#include <iomanip>

namespace NSPbm {
    using namespace NSPpredNA;
    using namespace std;

    map<int,char> baseTypeitoaMap = {{0,'A'}, {1, 'U'}, {2, 'G'}, {3, 'C'},
        {4,'a'}, {5, 't'}, {6, 'g'}, {7,'c'}, {-1, 'N'}};

    MotifGraph::MotifGraph(MotifAssigner* assigner, vector<int>& nodeIndexVec, double eneIn, const string& motifName) {
        fullGraph = assigner->nuGraph;
        name = motifName;
        nNode = nodeIndexVec.size();
        nFrag = 1;
        // sep.emplace_back(0); incorrect, the first motif node can not be guaranteed chain terminal.
        ene = eneIn;
        seq = new int[nNode];
        nodeIndex = new int[nNode];
        memcpy(nodeIndex, &nodeIndexVec[0], sizeof(int)*nNode);
        connectToDownstream = new bool[nNode];
        fixed = new bool[nNode];
        allNodes = new NuNode*[nNode];
        allEdges = new NuEdge*[nNode*nNode];
        geList.reserve(nNode*(nNode-1)/2);
        set<int> nodeIndexSet(nodeIndexVec.begin(), nodeIndexVec.end());
        pHelixEdge = assigner->pHelixEdge;
        pHelixNode2EdgeMap = assigner->pHelixNode2EdgeMap;  // Mapping from full-structure index to helix Edge
        pHelixNodeVec = assigner->pHelixNodeVec;  // vector of helix nodes
        pNode2HelixMap = assigner->pNode2HelixMap;  // map from nodeIndex to Helix numbering
        
        #ifdef DEBUG
            cout<<"[Info]MotifGraph: Calculating sepVec"<<endl;
        #endif
        sepVec.emplace_back(nodeIndex[0] ? getSep(nodeIndex[0],0) : 0);
        for(int i=0; i<nNode; i++) {
            allNodes[i] = fullGraph->allNodes[nodeIndex[i]];
            seq[i] = fullGraph->seq[nodeIndex[i]];
            if(i == nNode-1) {
                connectToDownstream[i] = false;
                // 3' end of the last fragment
                if(nodeIndex[nNode-1] == fullGraph->seqLen-1) sepVec.emplace_back(0);
                else sepVec.emplace_back(getSep(nodeIndex[nNode-1],fullGraph->seqLen-1));
            } else if(nodeIndexSet.count(nodeIndex[i]+1)) {
                connectToDownstream[i] = fullGraph->connectToDownstream[nodeIndex[i]];
                if(!connectToDownstream[i]) {  // chain break at i
                    sepVec.emplace_back(0);
                    nFrag+=1;
                    // check 5' break of the next fragment at nodeIndex[i]+1 coz nodeIndex[i+1] == nodeIndex[i]+1,
                    // 肯定是异链且是自然中断
                    sepVec.emplace_back(0);
                }
            } else {
                connectToDownstream[i] = false;  // the original downstream node is not included in subgraph
                // next motif node is nodeIndex[i+1], next full seq node is nodeIndex[i]+1
                // check 3' break
                sepVec.emplace_back(getSep(nodeIndex[i],nodeIndex[i+1]));
                // check 5' break of the next fragment
                nFrag+=1;
                sepVec.emplace_back(getSep(nodeIndex[i+1], nodeIndex[i]));
            }
            fixed[i] = fullGraph->fixed[nodeIndex[i]];
            for(int j=0; j<nNode; j++) {
                allEdges[i*nNode+j] = fullGraph->allEdges[nodeIndex[i]*fullGraph->seqLen+nodeIndex[j]];
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

    int MotifGraph::writeSeqFrag(ostream& outs) {
        string sf = name + ",";  // string of name + fragmented sequence
        for(int i=0; i<nNode-1; i++) {
            sf += baseTypeitoaMap[seq[i]];
            if(!connectToDownstream[i]) {
                sf+=",";
            }
        }
        sf+=baseTypeitoaMap[seq[nNode]];
        sf+="\n";
        outs << sf;
        return EXIT_SUCCESS;
    }

    int MotifGraph::writeNodeIndexInFull(ostream& outs) {
        string ni = name + ",";
        for(int i=0; i<nNode-1; i++) {
            ni += to_string(nodeIndex[i]) + ",";
        }
        ni += to_string(nodeIndex[nNode-1]) + "\n";
        outs << ni;
        return EXIT_SUCCESS;
    }

    int MotifGraph::writeMotif(ostream& outs, char label) {
        outs << left << setw(3) << name << " ";
        int iFull = 0;
        int i=0;
        while(i < nNode) {
            while(iFull < nodeIndex[i]) {
                outs << ".";
                iFull++;
            }
            outs <<  label ? label : baseTypeitoaMap[seq[nodeIndex[i]]];
            i++;
            iFull++;
        }
        while(iFull < fullGraph->seqLen) {
            outs << ".";
            iFull ++;
        }
        outs << "\n";
        return EXIT_SUCCESS;
    }


    int MotifGraph::writeSep(const string& outputFile) {
        if(nFrag != sepVec.size()/2) {
            throw "[Error] Corrupted sepVec";
        }
        ofstream ofs;
        ofs.open(outputFile, ios::out);
        if(!ofs.is_open()) {
            throw "[Error] writeSep: can not open file " + outputFile;
        }
        ofs << "fragment,3',5'\n";
        for(int i=0; i<nFrag; i++) {
            ofs << i <<","<<sepVec[i*2] << "," << sepVec[i*2+1] << '\n';
        }
        ofs.close();
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


    int MotifGraph::getSep(int ibeg, int iend) {
        int sep;
        if(ibeg < iend) {
            // 3' end of fragment, search to 5' end of the next fragment
            if(pNode2HelixMap->count(ibeg) &&
               pNode2HelixMap->count(ibeg+1)) {
                    // helix break: search the end of this helix at 5' side
                    sep = -1;  // sep = 0 is always false in the following if statement within the while-loop body
                    while(sep > -MOTIFGRAPH_SEP_SEARCH_MAX) {
                        if(!(pNode2HelixMap->count(ibeg+sep) && pNode2HelixMap->at(ibeg) == pNode2HelixMap->at(ibeg+sep))) {
                            return -(pNode2HelixMap->at(ibeg)*100+50+sep);  // within helix sep to helix 5' term 
                        }
                        sep--;
                    }
                    return -(pNode2HelixMap->at(ibeg)*100+1); // within helix exceeded sep to helix 5' term
            }
            for(int sep=0; sep<MOTIFGRAPH_SEP_SEARCH_MAX; sep++) {
                if(!fullGraph->connectToDownstream[ibeg+sep]) {
                    return sep;  // sep to chain break
                } else if(sep == iend-ibeg) {
                    return 300 + sep;  // sep to next fragment on the same chain
                } else if(pNode2HelixMap->count(ibeg+sep)) {
                    return -(pNode2HelixMap->at(ibeg+sep)*100+50+sep);  // sep to helix 5' term
                }
            }
            return sep;  // if no helix, chain break or the next fragment found, treat as chain break
        } else if(ibeg > iend) {
            // 5' end of fragment, search to 3' end of the last fragment
            if(pNode2HelixMap->count(ibeg) &&
               pNode2HelixMap->count(ibeg-1)) { //ibeg at helix break, search the end of this helix at 3' side
                sep = 1;    // sep = 0 is always false in the following if statement within the while-loop body
                while(sep < MOTIFGRAPH_SEP_SEARCH_MAX) {
                    if(!(pNode2HelixMap->count(ibeg+sep) && pNode2HelixMap->at(ibeg) == pNode2HelixMap->at(ibeg+sep))) {
                        return -(pNode2HelixMap->at(ibeg)*100+50+sep);  // within helix sep to helix 3' term 
                    }
                    sep++;
                }
                return -(pNode2HelixMap->at(ibeg)*100+99); // within helix exceeded sep to helix 3' term
            }
            for(int sep=0; sep>-MOTIFGRAPH_SEP_SEARCH_MAX; sep--) {
                if(!fullGraph->connectToDownstream[ibeg+sep-1]) {
                    return sep;  // sep to chain break
                } else if(sep == iend-ibeg) {
                    return 300 + sep;  // sep to last fragment on the same chain
                } else if(pNode2HelixMap->count(ibeg+sep)) {
                    return -(pNode2HelixMap->at(ibeg+sep)*100+50+sep);  // sep to helix 3' term
                }
            }
            return sep;  // if no helix, chain break or the last fragment found, treat as chain break
        } else {  // ibeg == iend
            throw "[Error] getSep: ibeg is identical to iend, ibeg=" + to_string(ibeg);
        }
    }


    MotifAssigner::MotifAssigner(NuGraph* nuGraphIn) {
        nuGraph = nuGraphIn;
        int seqLen = nuGraph->seqLen;
        baseStackingMap = getStackingSetMap();
        pHelixEdge = new set<NuEdge*>;
        for(int i=0; i<seqLen; i++) {
            for(int j=i+2; j<seqLen; j++) {
                NuEdge* curEdge = nuGraph->allEdges[i*seqLen+j];
                NuEdge* revEdge = nuGraph->allEdges[j*seqLen+i];
                if(isHelixByWC(i,j)) {
                    pHelixEdge->emplace(curEdge);
                    pHelixEdge->emplace(revEdge);
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
        pHelixNode2EdgeMap = new map<int, NuEdge*>;
        pHelixNodeVec = new vector<int>;
        pNode2HelixMap = new map<int,int>;
        #ifdef DEBUG
            cout<<"[Info]MotifAssigner: generating pHelixNode2EdgeMap"<<endl;
        #endif
        for(auto& it: *pHelixEdge) {
            if(pHelixNode2EdgeMap->count(it->indexA)) {
                if(it->indexA == pHelixNode2EdgeMap->at(it->indexA)->indexB &&
                   it->indexB == pHelixNode2EdgeMap->at(it->indexA)->indexA) {
                    continue; // skip reverse edge
                }
                string errInfo = "[ERROR] Duplicated Helix Node " + to_string(it->indexA) + " in edge" +\
                    to_string(it->indexA) + "-" + to_string(it->indexB) + "with" + \
                    to_string(pHelixNode2EdgeMap->at(it->indexA)->indexA) + "-" + \
                    to_string(pHelixNode2EdgeMap->at(it->indexA)->indexB) + "\n";
                cerr << errInfo;
                continue;
            }
            if(pHelixNode2EdgeMap->count(it->indexB)) {
                if(it->indexA == pHelixNode2EdgeMap->at(it->indexB)->indexB &&
                   it->indexB == pHelixNode2EdgeMap->at(it->indexB)->indexA) {
                    continue;  // skip reverse edge
                }
                string errInfo = "[ERROR] Duplicated Helix Node " + to_string(it->indexB) + " in edge" +\
                    to_string(it->indexA) + "-" + to_string(it->indexB) + "with" + \
                    to_string(pHelixNode2EdgeMap->at(it->indexB)->indexA) + "-" + \
                    to_string(pHelixNode2EdgeMap->at(it->indexB)->indexB) + "\n";
                cerr << errInfo;
                continue;
            }
            pHelixNode2EdgeMap->emplace(it->indexA, it);
            pHelixNode2EdgeMap->emplace(it->indexB, it);
            pHelixNodeVec->emplace_back(it->indexA);
            pHelixNodeVec->emplace_back(it->indexB);
        }
        sort(pHelixNodeVec->begin(),pHelixNodeVec->end());
        #ifdef DEBUG
            cout << "pHelixEdge has " << pHelixEdge->size() << " members" << endl;
            cout << "pHelixNodeVec has " << pHelixNodeVec->size() << " members" << endl;
            cout << "pHelixNode2EdgeMap has " << pHelixNode2EdgeMap->size() << " members" << endl;
            cout << "pHelixNodeVec" << endl;
            for(auto& it : *pHelixNodeVec) {
                cout << it << endl;
            }
        #endif
        int iHelix=1;
        int nHelixNode = pHelixNodeVec->size();
        #ifdef DEBUG
            cout<<"[Info]MotifAssigner: numbering helices"<<endl;
        #endif
        for(int i=0;i<nHelixNode-1;i++) {
            // check if current node has been treated
            // add nodes only when iterating the 5' half of helix
            if(pNode2HelixMap->count(pHelixNodeVec->at(i))) continue;  // jump over 3' half
            // add current node-iHelix mapping to node2HelixMap
            pNode2HelixMap->emplace(pHelixNodeVec->at(i), iHelix);
            NuEdge* curEdge = pHelixNode2EdgeMap->at(pHelixNodeVec->at(i));
            int pairNode = pHelixNodeVec->at(i) == curEdge->indexA ? curEdge->indexB : curEdge->indexA;
            pNode2HelixMap->emplace(pairNode, iHelix);
            
            // check if helix continues
            NuEdge* nextEdge = pHelixNode2EdgeMap->at(pHelixNodeVec->at(i+1));
            int pairNodeNext = pHelixNodeVec->at(i+1) == nextEdge->indexA ? nextEdge->indexB : nextEdge->indexA;
            #ifdef DEBUG
                if(pHelixNodeVec->at(i) == 0) {
                    cout << "node 0 helix next node " << pHelixNodeVec->at(i+1) << ", next pairNode " << pairNodeNext << endl;
                    cout << "full graph connectTodownStream[0]=" << boolalpha << nuGraph->connectToDownstream[0] <<
                          ", full graph connectToDownStream[pairNodeNext]=" << nuGraph->connectToDownstream[pairNodeNext] << noboolalpha << endl;
                }
            # endif
            if(pHelixNodeVec->at(i+1)-pHelixNodeVec->at(i) == 1 && nuGraph->connectToDownstream[pHelixNodeVec->at(i)] &&
               pairNode-pairNodeNext == 1 && nuGraph->connectToDownstream[pairNodeNext]) {
                continue;
            } else {
                iHelix+=1;
            }
        }
        #ifdef DEBUG
            // print ihelix nodes
            cout << "pNode2HelixMap:" <<endl;
            for(auto& it : *pNode2HelixMap) {
                cout << "node: " << it.first << ", helix: " << it.second << endl;
            }            
        #endif
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
     * a stronger stacking&basepair weight than edge ij.
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
            for(auto& k : *(baseStackingMap->at(i))) {
                if(min(nuGraph->allEdges[k*seqLen+j]->weight,
                       nuGraph->allEdges[j*seqLen+k]->weight) < curWeight &&
                   min(nuGraph->allEdges[k*seqLen+i]->weight,
+                      nuGraph->allEdges[i*seqLen+k]->weight) < curWeight) {
                    return true;
                }
            }
        }
        if(baseStackingMap->count(j)) {
            for(auto& k : *(baseStackingMap->at(j))) {
                if(min(nuGraph->allEdges[k*seqLen+i]->weight,
                       nuGraph->allEdges[i*seqLen+k]->weight) < curWeight &&
                   min(nuGraph->allEdges[k*seqLen+j]->weight,
+                      nuGraph->allEdges[j*seqLen+k]->weight) < curWeight) {
                    return true;
                }
            }
        }
        return false;
    }


    /**
     * @brief check if all nodes in **nodeList** belong to the
     * same helix denoted by **iHelix**
     * 
     * @param i 
     * @param iHelix 
     * @param nodeList 
     * @return true 
     * @return false 
     */
    bool MotifAssigner::isAllNodeInHelix(int iHelix, set<int>* nodeList) {
        for(const int& it: *nodeList) {
            try {
                if(pNode2HelixMap->at(it) != iHelix) {
                    return false;
                }
            } catch (out_of_range) {
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


    void MotifAssigner::bySeed(string&& outPath, string&& outPDBprefix) {
        vector<NuEdge*> seedEdge;  // take non-neighboring nonWC pairs as seed
        int seqLen = nuGraph->seqLen;
        for(int i=0; i<seqLen; i++) {
            for(int j=i+2; j<seqLen; j++) {  // Exclude sequence neighbor
                NuEdge* curEdge = nuGraph->allEdges[i*seqLen+j];
                NuEdge* revEdge = nuGraph->allEdges[j*seqLen+i];

                // Exclude helix WC Edge in seed search
                if(curEdge->isWC() && pHelixEdge->count(curEdge)) {
                    continue;
                }

                // Exclude implicit edge
                if(isImplicitEdge(i,j)) {
                    continue;
                }

                // Exclude edge between nodes of a same helix
                if(pNode2HelixMap->count(i) && pNode2HelixMap->count(j) &&
                   pNode2HelixMap->at(i) == pNode2HelixMap->at(j)) {
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

        // set<int> WCNodeSet;  // [deprecated] Set of nodes that have only WC non-neighbor basepair
        // for(auto& it: nnbNeighborMap) {
        //     if(it.second->size() == 1 && 
        //     nuGraph->allEdges[it.first*seqLen + *(it.second->begin())]->isWC()) {
        //         WCNodeSet.emplace(it.first);
        //     }
        // }

        int nSeed = seedEdge.size();
        ofstream outCSV;
        string csvFileName = outPath + "/" + outPDBprefix + "-Seeds" + ".csv";
        outCSV.open(csvFileName, ios::out);
        if(!outCSV.is_open()) {
            throw "[Error] Fail to open file" + csvFileName;
        }
        outCSV << "indexA,indexB,weight,isWC" << endl;
        for(int i=0; i<nSeed; i++) {
            set<int> selNode{seedEdge[i]->indexA, seedEdge[i]->indexB};
            cout << "grow from seed node: " << seedEdge[i]->indexA << "," << seedEdge[i]->indexB << endl;
            outCSV << seedEdge[i]->indexA << ","
                   << seedEdge[i]->indexB << ","
                   << seedEdge[i]->weight << ","
                   << seedEdge[i]->isWC() << endl;
            set<int>* pCurConnect = new set<int>;
            set<int>* pLastConnect = new set<int>, *ptmp;
            set<NuEdge*> selEdge = {seedEdge[i]};
            set<int> helixTermNodes = set<int>{};
            // map<int,int>* nodeGrowCt = new map<int,int>;
            // set<int>* pToExclude = new set<int>;
            // set<int>* pPostExclude = new set<int>;
            #ifdef DEBUG
                cout << "Initialize pCurConnect set" << endl;
            #endif
            for(auto& it:selNode) {
                #ifdef DEBUG
                    cout << "[Info]BySeed: init round, base " << it << endl;
                #endif
                // connect all nnb neighbors (non-sequential-neighboring interactions)
                // set_union(pCurConnect->begin(), pCurConnect->end(), 
                //           nnbNeighborMap[it]->begin(), nnbNeighborMap[it]->end(),
                //           inserter(*pCurConnect, pCurConnect->begin()));
                bool bCurPureHelix = false;  // flags for current node in pure helix
                    // current node is thought pure helix if all nnb edges are in the same helix
                bool bGrowPureHelix = false;   // flags for grow node in the same pure helix
                    // grow node is thought pure helix if all nb and nnb edges are in the same helix as the current node
                int curHelix = -1;
                if(pNode2HelixMap->count(it)) {
                    curHelix = pNode2HelixMap->at(it);
                    NuEdge* curHelixEdge = pHelixNode2EdgeMap->at(it);
                    int pairNode = curHelixEdge->indexA == it? curHelixEdge->indexB: curHelixEdge->indexA;
                    bCurPureHelix = isAllNodeInHelix(curHelix, nnbNeighborMap[it]) && isAllNodeInHelix(curHelix, nnbNeighborMap[pairNode]);
                    #ifdef DEBUG
                        cout << "current pair: " << it << "-" << pairNode <<endl;
                        cout << "curHelix: " << curHelix << endl;
                        cout << "Is current helix pure: " << bCurPureHelix << endl;
                    #endif
                }
                for(const int& it1: *(nnbNeighborMap[it])) {
                    bool bWCEdge = nuGraph->allEdges[it*seqLen + it1]->isWC();  // allow WC edge
                    if(!bWCEdge && bCurPureHelix && pNode2HelixMap->count(it1) && pNode2HelixMap->at(it1) == curHelix) {
                        bGrowPureHelix = isAllNodeInHelix(curHelix, nbNeighborMap[it1]) && \
                                         isAllNodeInHelix(curHelix, nnbNeighborMap[it1]);
                        // include nodes but mark as helix terminal when all nnb edges of **it** is in curhelix and all edges of **it1** is in curhelix
                        if(bGrowPureHelix) helixTermNodes.emplace(it1);
                    }
                    pCurConnect->emplace(it1);
                }
                // connect nb neighbors, set boundary on helix
                for(const int& it2: *(nbNeighborMap[it])) {
                    // [deprecated] if current node is WC-only node and the nb neighbor is also WC-only node, then grow will stop,
                    // so the nb neighbor excluded
                    // if(WCNodeSet.count(it2) && WCNodeSet.count(it)) {
                    //     continue;
                    // }                    
                    
                    // nbEdge can not be WC, no need to check bWCEdge for it-it2
                    if(bCurPureHelix && pNode2HelixMap->count(it2) && pNode2HelixMap->at(it2) == curHelix) {
                        bGrowPureHelix = isAllNodeInHelix(curHelix, nbNeighborMap[it2]) && \
                                         isAllNodeInHelix(curHelix, nnbNeighborMap[it2]);
                        // include nodes but mark as helix terminal when all nnb edges of **it** is in curHelix and all edges of **it2** is in curHelix
                        if(bGrowPureHelix) helixTermNodes.emplace(it2);
                    }
                    pCurConnect->emplace(it2);
                }
                #ifdef DEBUG
                    cout<< "pCurConnect: ";
                    for(auto& it : *pCurConnect) {
                        cout << it <<",";
                    }
                    cout << endl;
                    cout<< "helixTermNodes: ";
                    for(auto& it : helixTermNodes) {
                        cout << it <<",";
                    }
                    cout << endl;
                #endif
            }
            // nodeGrowCt->clear();
            // pToExclude->clear();
            // pPostExclude->clear();
            #ifdef DEBUG
                cout << "iterate until pCurConnect nodes are already included in selNode" << endl;
                int nround=1;
            #endif
            while(!includes(selNode.begin(), selNode.end(), pCurConnect->begin(), pCurConnect->end())) {
                set_union(selNode.begin(), selNode.end(),
                          pCurConnect->begin(), pCurConnect->end(),
                          inserter(selNode, selNode.begin()));
                //TODO selEdge
                pLastConnect->clear();
                // initialize the connected nodes in next round according to pCurConnect, use pLastConnect temporarily
                for(auto& it:*pCurConnect) {
                    #ifdef DEBUG
                        cout << "[Info]BySeed: round " << nround <<", base " << it << endl;
                    #endif
                    if(helixTermNodes.count(it)) { // do not grow at helix terminal nodes, except WC edge
                        #ifdef DEBUG
                            cout << it << " is at helix terminal" << endl;
                        #endif
                        for(const int& it1: *(nnbNeighborMap[it])) {
                            bool bWCEdge = nuGraph->allEdges[it*seqLen + it1]->isWC();  // allow WC edge
                            if(bWCEdge) {
                                pLastConnect->emplace(it1);
                                helixTermNodes.emplace(it1);
                                break; // it is not possible to have 2 WC edge on a same node
                            }
                        }
                        continue;
                    }
                    bool bCurPureHelix = false;  // flags for current node in pure helix
                        // current node is thought pure helix if all nnb edges are in the same helix
                    bool bGrowPureHelix = false;   // flags for grow node in the same pure helix
                        // grow node is thought pure helix if all nb and nnb edges are in the same helix as the current node
                    int curHelix = -1;
                    if(pNode2HelixMap->count(it)) {
                        curHelix = pNode2HelixMap->at(it);
                        NuEdge* curHelixEdge = pHelixNode2EdgeMap->at(it);
                        int pairNode = curHelixEdge->indexA == it? curHelixEdge->indexB: curHelixEdge->indexA;
                        bCurPureHelix = isAllNodeInHelix(curHelix, nnbNeighborMap[it]) && isAllNodeInHelix(curHelix, nnbNeighborMap[pairNode]);
                        #ifdef DEBUG
                            cout << "current pair: " << it << "-" << pairNode <<endl;
                            cout << "curHelix: " << curHelix << endl;
                            cout << "Is current helix pure: " << bCurPureHelix << endl;
                        #endif
                    }
                    for(const int& it1: *(nnbNeighborMap[it])) {
                        bool bWCEdge = nuGraph->allEdges[it*seqLen + it1]->isWC();  // allow WC edge
                        if(!bWCEdge && bCurPureHelix && pNode2HelixMap->count(it1) && pNode2HelixMap->at(it1) == curHelix) {
                            bGrowPureHelix = isAllNodeInHelix(curHelix, nbNeighborMap[it1]) && \
                                             isAllNodeInHelix(curHelix, nnbNeighborMap[it1]);
                            // exclude nodes when all nnb edges of **it** is in curhelix and all edges of **it1** is in curhelix
                            if(bGrowPureHelix) helixTermNodes.emplace(it1);
                        }
                        pLastConnect->emplace(it1);
                    }
                    // set_union(pLastConnect->begin(), pLastConnect->end(), 
                    //           nnbNeighborMap[it]->begin(), nnbNeighborMap[it]->end(),
                    //           inserter(*pLastConnect, pLastConnect->begin()));
                    for(const int& it2: *(nbNeighborMap[it])) {
                        // if(WCNodeSet.count(it2) && WCNodeSet.count(it)) { //[deprecated]
                        //     continue;
                        // }
                        if(bCurPureHelix && pNode2HelixMap->count(it2) && pNode2HelixMap->at(it2) == curHelix) {
                            bGrowPureHelix = isAllNodeInHelix(curHelix, nbNeighborMap[it2]) && \
                                             isAllNodeInHelix(curHelix, nnbNeighborMap[it2]);
                            // exclude nodes when all nnb edges of **it** is in curHelix and all edges of **it2** is in curHelix
                            if(bGrowPureHelix) helixTermNodes.emplace(it2);
                        }
                        pLastConnect->emplace(it2);
                    }

                    // exclude nodes when all nnb edges of **it** is in helix1 and all edges of **it2** is in helix1
                    // if(pNode2HelixMap->count(it)) {
                    //     int curHelix = pNode2HelixMap->at(it);
                    //     bool itPureHelix = true;
                    //     for(const int& it3: *(nnbNeighborMap[it])) {
                    //         try {
                    //             if(pNode2HelixMap->at(it3) != curHelix) {
                    //                 itPureHelix = false;
                    //                 break;
                    //             }
                    //         } catch (out_of_range) {
                    //             itPureHelix = false;
                    //             break;
                    //         }
                    //     }
                    //     if(!itPureHelix) continue;

                    //     set<int>* pExclude = new set<int>;
                    //     for(const int& it4 : *pLastConnect) {
                    //         try {
                    //             if(pNode2HelixMap->at(it4) != curHelix) {
                    //                 continue;
                    //             }
                    //         } catch (out_of_range) {
                    //             continue;
                    //         }
                    //         bool it4PureHelix = true;
                    //         for(const int& it5: *(nnbNeighborMap[it4])) {
                    //             try {
                    //                 if(pNode2HelixMap->at(it5) != curHelix) {
                    //                     it4PureHelix = false;
                    //                     break;
                    //                 }
                    //             } catch (out_of_range) {
                    //                 it4PureHelix = false;
                    //                 break;
                    //             }
                    //         }
                    //         for(const int& it6: *(nbNeighborMap[it4])) {
                    //             try {
                    //                 if(pNode2HelixMap->at(it6) != curHelix) {
                    //                     it4PureHelix = false;
                    //                     break;
                    //                 }
                    //             } catch (out_of_range) {
                    //                 it4PureHelix = false;
                    //                 break;
                    //             }
                    //         }
                    //         if(it4PureHelix) pExclude->emplace(it4);
                    //     }

                    //     set<int>* pResult = new set<int>, *pTmp;
                    //     set_difference(pLastConnect->begin(), pLastConnect->end(),
                    //                    pExclude->begin(), pExclude->end(),
                    //                    inserter(*pResult, pResult->begin()));
                    //     pTmp = pLastConnect;
                    //     pLastConnect = pResult;
                    //     delete pTmp; 
                    //     delete pExclude;                   
                    // }
                    #ifdef DEBUG
                        cout<< "current Connect: ";
                        for(auto& it : *pLastConnect) {
                            cout << it <<",";
                        }
                        cout << endl;
                        cout<< "helixTermNodes: ";
                        for(auto& it : helixTermNodes) {
                            cout << it <<",";
                        }
                        cout << endl;
                    #endif
                }
                // point pLastConnect to the last round pCurConnect, and pCurConnect to the next round connected nodes
                ptmp = pCurConnect;
                pCurConnect = pLastConnect;
                pLastConnect = ptmp;
                #ifdef DEBUG
                    nround++;
                #endif
            }
            #ifdef DEBUG
                cout << "[info]bySeed: clearing allocated memory" << endl;
            #endif
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

            MotifGraph* pMotif = new MotifGraph(this, nodeIndexVec, ene);

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
        outCSV.close();

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

    int MotifAssigner::writeSeqFrag(ostream& outs) {
        outs << "Motif,seqFrag" << endl;
        for(auto& it: motifs) {
            it->writeSeqFrag(outs);
        }
        return EXIT_SUCCESS;
    }

    int MotifAssigner::writeNodeIndexInFull(ostream& outs) {
        outs << "Motif, Index" << endl;
        for(auto& it: motifs) {
            it->writeNodeIndexInFull(outs);
        }
        return EXIT_SUCCESS;
     }

    int MotifAssigner::writeMotif(ostream& outs, string& ssSeq) {
        string seqStr, connectStr;
        for(int i=0; i<nuGraph->seqLen; i++) {
            seqStr += baseTypeitoaMap[nuGraph->seq[i]];
            connectStr += nuGraph->connectToDownstream[i]? "-":"|";
        }
        outs << "seq " << seqStr << endl;
        outs << "sec " << ssSeq << endl;
        outs << "cnt " << connectStr << endl;
        int nMotif = motifs.size();
        for(int i=0; i<nMotif; i++) {
            string name0 = motifs[i]->name;
            motifs[i]->setName("m" + to_string(i));
            motifs[i]->writeMotif(outs, 'A'+i%26);
            motifs[i]->setName(name0);
        }
        return EXIT_SUCCESS;
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
        delete pHelixEdge;
        delete pHelixNode2EdgeMap;
        delete pHelixNodeVec;
        delete pNode2HelixMap;
    }
}