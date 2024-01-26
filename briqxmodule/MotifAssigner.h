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

#define MOTIFASSIGNER_SEED_CUTOFF 0.5
#define MOTIFASSIGNER_GROW_CUTOFF 0.5
#define MOTIFASSIGNER_RECALL_CUTOFF 0.25

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
            virtual ~MotifGraph() {
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
            vector<MotifGraph*> motifs;
            MotifAssigner(NuGraph* nuGraphIn) {
                nuGraph = nuGraphIn;
            };
            void bySeed();  // assign motif by greedy algorithm starting from seed BPs
            void byLouvain();

            virtual ~MotifAssigner();
    };
}
#endif /* BRIQXMODULE_MOTIFASSIGNER_H_ */