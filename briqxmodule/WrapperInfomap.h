/**
 * @file WrapperInfomap.h
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief Wrap Infomap to be used with NuGraph
 * @version 0.1
 * @date 2023-12-28
 * 
 * @copyright Copyright (c) 2023 XLAB
 * 
 */

#ifndef BRIQXMODULE_WRAPPERINFOMAP_H_
#define BRIQXMODULE_WRAPPERINFOMAP_H_

#include "predNA/NuGraph.h"
// #include "io/Network.h"

namespace NSPbm {
    using namespace NSPpredNA;
    using namespace std;

    class WrapperInfomap {
        public:
            NuGraph* nuGraph;
            // infomap::Network infoGraph;

            /**
             * @brief Construct a new Wrapper Infomap object with given NuGraph
             * 
             * @param nuGraphIn
             */
            WrapperInfomap(NuGraph& nuGraphIn);

            void run();

            void print() const;

            /**
             * @brief Write NuGraph in link list format
             * 
             */
            void writeLinkList(ostream &out) const;

            /**
             * @brief Read infomap generated community info and mapping to NuGraph/create new NuGraphs
             * 
             */
            void readCommunities();
    };
}
#endif /* BRIQXMODULE_WRAPPERINFOMAP_H_ */
