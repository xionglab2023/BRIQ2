/**
 * @file testMoveIO.cpp
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief test IO performance of NuMoveSet with txt and binary data sources
 * @version 0.1
 * @date 2024-01-09
 * 
 * @copyright Copyright (c) 2024 XLAB
 * 
 */

#include "predNA/NuMoveSet.h"
#include <sys/time.h>

using namespace NSPdataio;
using namespace NSPpredNA;
using namespace std;

int main() {
    struct timeval start, end;
    double timeTxt, timeBinary;
    cout << "Creating NuMoveSet from txt" << endl;
    gettimeofday(&start, NULL);  // gettimeofday 计现实时间，clock() 计cpu时钟数，受sleep() 和多线程影响
    auto * set1 = new NuPairMoveSetLibrary(false);
    gettimeofday(&end, NULL);
    timeTxt = end.tv_sec - start.tv_sec + (double)(end.tv_usec-start.tv_usec) /1e6;
    cout << "time consumed: " << timeTxt << " seconds" << endl;
    cout << "first int in each vector of nbMoveList[6][9] moveIndexList:" << endl;
    for(int i=0; i<20; i++) {
        cout << set1->nbMoveList[6][9]->moveIndexList[i][0] << " ";
    }
    cout << endl;
    cout << "first int in each vector of revNbMoveList[6][9] moveIndexList:" << endl;
    for(int i=0; i<20; i++) {
        cout << set1->revNbMoveList[6][9]->moveIndexList[i][0] << " ";
    }
    cout << endl;
    cout << "first int in each vector of nnbMoveList[6][9] moveIndexList:" << endl;
    for(int i=0; i<20; i++) {
        cout << set1->nnbMoveList[6][9]->moveIndexList[i][0] << " ";
    }
    cout << endl;
    delete set1;
    
    cout << "Creating NuMoveSet from binary" << endl;
    gettimeofday(&start, NULL);
    set1 = new NuPairMoveSetLibrary();
    gettimeofday(&end, NULL);
    timeBinary = end.tv_sec - start.tv_sec + (double)(end.tv_usec-start.tv_usec) /1e6;
    cout << "time consumed: " << timeBinary << " seconds" << endl;
    cout << "first int in each vector of nbMoveList[6][9] moveIndexList:" << endl;
    for(int i=0; i<20; i++) {
        cout << set1->nbMoveList[6][9]->moveIndexList[i][0] << " ";
    }
    cout << endl;
    cout << "first int in each vector of revNbMoveList[6][9] moveIndexList:" << endl;
    for(int i=0; i<20; i++) {
        cout << set1->revNbMoveList[6][9]->moveIndexList[i][0] << " ";
    }
    cout << endl;
    cout << "first int in each vector of nnbMoveList[6][9] moveIndexList:" << endl;
    for(int i=0; i<20; i++) {
        cout << set1->nnbMoveList[6][9]->moveIndexList[i][0] << " ";
    }
    cout << endl;
    cout << "Reading Binary is " << timeTxt/timeBinary << "x boosting than reading txt." << endl;
    delete set1;
}