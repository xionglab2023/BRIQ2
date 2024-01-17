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

 namespace NSPbm {
    using namespace NSPpredNA;
    using namespace std;

    class MotifAssigner {
        public:
            NuGraph* nuGraph;
            vector<NuGraph*> motifs;
            MotifAssigner(NuGraph* nuGraphIn);
            void bySeed();
            void byLouvain();

            virtual ~MotifAssigner();
    };
 }
#endif /* BRIQXMODULE_MOTIFASSIGNER_H_ */