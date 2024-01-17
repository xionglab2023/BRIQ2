/**
 * @file WrapperInfomap.cpp
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief Implemtations of WrapperInfomap
 * @version 0.1
 * @date 2023-12-28
 * 
 * @copyright Copyright (c) 2023 XLAB
 * 
 */

#include "briqxmodule/WrapperInfomap.h"

namespace NSPbm {
    WrapperInfomap::WrapperInfomap(NuGraph& nuGraphIn) {
        nuGraph = &nuGraphIn;
    }

    void WrapperInfomap::run() {
        cout << "To be Implemented" << endl;
        return;
    }

    void WrapperInfomap::print() const {
        cout << "To be Implemented" << endl;
        return;
    }

    void WrapperInfomap::writeLinkList(ostream &out) const {
        out << "# Network derived from BRIQX NuGraph" << endl;
        out << "# " << endl;
        out << "# source target weight" << endl;
        int le = nuGraph->seqLen * nuGraph->seqLen;
        for (int i=0; i< le; i++) {
            NuEdge* nuEdge = nuGraph->allEdges[i];
            if(nuEdge->weight != 0) {
                out << nuEdge->indexA <<" "<< nuEdge->indexB <<" "<< nuEdge->weight <<endl;
            }
        }
    }
}