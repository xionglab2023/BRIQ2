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

#define MOTIFASSIGNER_SEED_CUTOFF -8.5
#define MOTIFASSIGNER_GROW_CUTOFF -5.5
#define MOTIFASSIGNER_RECALL_CUTOFF -7.0
#define MOTIFASSIGNER_STACK_WEIGHT_MAX -1e-6
#define MOTIFASSIGNER_BP_WEIGHT_MAX -12
#define DEBUG

namespace NSPbm {
    using namespace NSPpredNA;
    using namespace std;

    /**
     * @brief Class to store a subgraph that represents assigned motif, it holds an array recording the
     * initial index of nodes in the original NuGraph.
     * 
     */
    class MotifGraph {
        public:
            int nNode;  // L
            double ene;
            NuGraph* fullGraph;
            int* seq;
            int* nodeIndex;  // fullGraph index of the NuNodes in MotifGraph
            bool* connectToDownstream;
            bool* fixed;  // for compatiblity with graphInfo
            NuNode** allNodes;  // L nodes
            NuEdge** allEdges;  // L*L edges
            vector<NuEdge*> geList;  //L*(L-1)/2 edges
            graphInfo* motifInfo;
            // map<array<NuNode*,2>,NuEdge*> node2EdgeMap;

            /**
             * @brief Construct a new Motif Graph object
             * 
             * @param[in] GraphIn The original NuGraph.
             * @param[in] NodeIndexVec Indexes of selected nodes in GraphIn allNodes (and seq)
             */
            MotifGraph(NuGraph* GraphIn, vector<int>& NodeIndexVec, double ene=0.0);
            int writePDB(const string& outputFile);
            
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
    };


    class MotifAssigner {
        public:
            NuGraph* nuGraph;
            set<NuEdge*> helixEdge;
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
            #ifdef DEBUG
            void bySeed(string&&, string&&);
            #else
            void bySeed();  // assign motif by greedy algorithm starting from seed BPs
            #endif
            void byLouvain();

            virtual ~MotifAssigner();
    };
}
#endif /* BRIQXMODULE_MOTIFASSIGNER_H_ */