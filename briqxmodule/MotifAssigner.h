/**
 * @file MotifAssigner.h
 * @author Klark Chen  (klarkchen@ustc.edu.cn)
 * @brief Worker class to assign motif from given NuGraph, with several different algorithms
 * @version undef
 * @date 2024-01-03
 * 
 * @copyright Copyright (c) 2024 XLAB
 * 
 */

#ifndef BRIQXMODULE_MOTIFASSIGNER_H_
#define BRIQXMODULE_MOTIFASSIGNER_H_
 
#include "predNA/NuGraph.h"

#define MOTIFASSIGNER_SEED_CUTOFF -3
#define MOTIFASSIGNER_GROW_CUTOFF -1.5
#define MOTIFASSIGNER_BP_WEIGHT_MAX -4.6
#define MOTIFGRAPH_SEP_SEARCH_MAX 10  // can not exceed 50
// #define MOTIFASSIGNER_RECALL_CUTOFF -7.0
// #define DEBUG

namespace NSPbm {
    using namespace NSPpredNA;
    using namespace std;

    class MotifAssigner;
    /**
     * @brief Class to store a subgraph that represents assigned motif, it holds an array recording the
     * initial index of nodes in the original NuGraph.
     * 
     */
    class MotifGraph {
        public:
            int nNode;  // L
            int nFrag;  // number of fragments
            double ene;
            NuGraph* fullGraph;
            int* seq;
            int* nodeIndex;  // fullGraph index of the NuNodes in MotifGraph
            bool* connectToDownstream;
            bool* fixed;  // for compatiblity with graphInfo
            NuNode** allNodes;  // L nodes
            NuEdge** allEdges;  // L*L edges
            set<NuEdge*>* pHelixEdge;
            map<int,NuEdge*>* pHelixNode2EdgeMap;  // Mapping from full-structure index to helix Edge
            vector<int>* pHelixNodeVec;
            map<int,int>* pNode2HelixMap;  // map from nodeIndex to Helix numbering
            string name;
            vector<NuEdge*> geList;  //L*(L-1)/2 edges
            graphInfo* motifInfo;
            vector<int> sepVec;  // nFrag*2 sized array, recording the number of bases to the last/next fragment.
            // map<array<NuNode*,2>,NuEdge*> node2EdgeMap;

            /**
             * @brief Construct a new Motif Graph object
             * 
             * @param[in] GraphIn The original NuGraph.
             * @param[in] NodeIndexVec Indexes of selected nodes in GraphIn allNodes (and seq), the indexes
             *      may not be continuous but must be sorted.
             */
            MotifGraph(MotifAssigner* assigner, vector<int>& NodeIndexVec, double ene=0.0,
                 const string& motifName="defaultMotif");
            
            void setName(const string& motifName) {
                name = motifName;
            }
            int writePDB(const string& outputFile);
            int writeSeqFrag(ostream& outs);  // write sequence with fragments seperated by comma
            int writeNodeIndexInFull(ostream& outs);  // write indices in full structure for nodes included in motif
            int writeMotif(ostream& outs, char label=0);  // write motif info in briqxmotif format
            /**
             * @brief write sepVec in CSV format. It writes one file per MotifGraph
             * 
             * @param outputFile 
             * @return int 
             */
            int writeSep(const string& outputFile);
            
            /**
             * @brief Test if **this** is identical to **other**, where identical means same **fullGraph**,
             * same **nNode** and **same** nodeIndex
             * 
             * @param other 
             * @return true **this** is identical to **other**
             * @return false 
             */
            bool operator==(const MotifGraph& other);

            /**
             * @brief Test if **other** is included in **this**
             * 
             * @param other 
             * @return true 
             * @return false 
             */
            bool includes(const MotifGraph& other);
            virtual ~MotifGraph() {
                delete[] seq;
                delete[] nodeIndex;
                delete[] connectToDownstream;
                delete[] allNodes;
                delete[] allEdges;
                delete motifInfo;
            }
        private:
            /**
             * @brief Determine sep from index ibeg to iend according to sep code rules. Priority ranking of
             * sep types: helix, chain break, motif break.
             * 
             * @param ibeg: the full-structure index of the begining node of search (the break point)
             * @param iend: the full-structure index of the end node of search  (the next break point)
             * @return int: sep value 
             */
            int getSep(int ibeg, int iend);
    };


    class MotifAssigner {
        public:
            NuGraph* nuGraph;
            set<NuEdge*>* pHelixEdge;
            map<int,NuEdge*>* pHelixNode2EdgeMap;  // Mapping from full-structure index to helix Edge
            vector<int>* pHelixNodeVec;  // vector of helix nodes
            map<int,int>* pNode2HelixMap;  // map from nodeIndex to Helix numbering
            map<int,set<int>*>* baseStackingMap;
            vector<MotifGraph*> motifs;
            MotifAssigner(NuGraph* nuGraphIn);

            /**
             * @brief Write edge weights in csv format, with 5 columns:
             * nodeIndex1, nodeIndex2, weight, isWC, isNB
             * 
             * @param outCSV the output stream object, must in txt mode
             * @return int execute state code
             */
            int writeEdgeWeight(ostream& outCSV);  // write edge weights in csv format, with feature columns
            int writeSeqFrag(ostream& out);  // write motif sequeces with fragmentation
            int writeNodeIndexInFull(ostream& outs);  // write indices in full structure for nodes included in motifs
            int writeMotif(ostream& outs, string& ssSeq);  // write bm file (briqxMotif)
            void bySeed(string&&, string&&);
            void byLouvain();

            virtual ~MotifAssigner();
        private:
            bool isHelixByWC(int i, int j);
            bool isHelixByStacking(int i, int j);
            map<int,set<int>*>* getStackingSetMap();  // returns a map on the stack, requires manual release;
            bool isImplicitEdge(int i, int j);
            bool isAllNodeInHelix(int iHelix, set<int>* nodeList);
    };
}
#endif /* BRIQXMODULE_MOTIFASSIGNER_H_ */