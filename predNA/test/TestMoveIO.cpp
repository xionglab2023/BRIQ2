/**
 * @file testMoveIO.cpp
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief test IO performance of NuMoveSet with txt and binary data soutces
 * @version 0.1
 * @date 2024-01-09
 * 
 * @copyright Copyright (c) 2024 XLAB
 * 
 */

#include "predNA/NuMoveSet.h"
#include <time.h>

using namespace NSPdataio;
using namespace NSPpredNA;
using namespace std;

int main() {
    clock_t start, end;
    double timeTxt, timeBinary;
    cout << "Creating NuMoveSet from txt" << endl;
    start = clock();
    auto * set1 = new NuPairMoveSetLibrary(false);
    end = clock();
    timeTxt = (double)(end-start) / CLOCKS_PER_SEC;
    cout << "time consumed: " << timeTxt << " seconds" << endl;
    delete set1;
    
    cout << "Creating NuMoveSet from binary" << endl;
    start = clock();
    set1 = new NuPairMoveSetLibrary();
    end = clock();
    timeBinary = (double)(end-start) / CLOCKS_PER_SEC;
    cout << "time consumed: " << timeBinary << " seconds" << endl;
    cout << "Reading Binary is " << timeTxt/timeBinary << "x boosting than reading txt." << endl;
    delete set1;
}