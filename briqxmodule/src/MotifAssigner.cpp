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

namespace NSPbm {
    using namespace NSPpredNA;
    using namespace std;

    MotifAssigner::MotifAssigner(NuGraph* nuGraphIn) {
        nuGraph = nuGraphIn;
    }

    void MotifAssigner::bySeed() {
        return;
    }

    void MotifAssigner::byLouvain() {
        return;
    }

    MotifAssigner::~MotifAssigner() {
        //TODO
    }
}