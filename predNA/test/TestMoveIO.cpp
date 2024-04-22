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

int main(int argc, char** argv) {


    string libType = string(argv[1]);

    struct timeval start, end;
    double timeTxt, timeBinaryTable, timeBinaryCache;
    cout << "Creating NuMoveSet from txt" << endl;
    gettimeofday(&start, NULL);  // gettimeofday 计现实时间，clock() 计cpu时钟数，受sleep() 和多线程影响
    auto * set1 = new NuPairMoveSetLibrary(libType, false);
    gettimeofday(&end, NULL);
    timeTxt = end.tv_sec - start.tv_sec + (double)(end.tv_usec-start.tv_usec) /1e6;
    cout << "time consumed: " << timeTxt << " seconds" << endl;

    set1->printMoveLibInfo();

    cout << "first int in each vector of nbMoveList[6][9] moveIndexList:" << endl;
    for(int i=0; i<20; i++) {
        cout << set1->getMoveSet(6,9,1)->moveIndexList[i][0] << " ";
    }
    cout << endl;
    cout << "first int in each vector of revNbMoveList[6][9] moveIndexList:" << endl;
    for(int i=0; i<20; i++) {
        cout << set1->getMoveSet(6,9,-1)->moveIndexList[i][0] << " ";
    }
    cout << endl;
    cout << "first int in each vector of nnbMoveList[6][9] moveIndexList:" << endl;
    for(int i=0; i<20; i++) {
        cout << set1->getMoveSet(6,9,2)->moveIndexList[i][0] << " ";
    }
    cout << endl;
    cout << "oi sphereKeyMap1000 first 5 key-val pairs: " << endl;
    int ii = 0;
		for(auto& it:set1->oi->getSKM1000()){
			cout << it.first << ":" << it.second << endl;
			ii++;
            if(ii>4) break;
		}
    set1->dump();
    delete set1;
    
    /*
    cout << "Creating NuMoveSet from binaryTable" << endl;
    gettimeofday(&start, NULL);
    set1 = new NuPairMoveSetLibrary(true,2);
    gettimeofday(&end, NULL);
    timeBinaryTable = end.tv_sec - start.tv_sec + (double)(end.tv_usec-start.tv_usec) /1e6;
    cout << "time consumed: " << timeBinaryTable << " seconds" << endl;
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
    cout << "oi sphereKeyMap1000 first 5 key-val pairs: " << endl;
    ii = 0;
	for(auto& it:set1->oi->getSKM1000()){
		cout << it.first << ":" << it.second << endl;
		ii++;
           if(ii>4) break;
	}
    cout << "Reading BinaryTable is " << timeTxt/timeBinaryTable << "x boosting than reading txt." << endl;
    delete set1;
    */

    cout << "Creating NuMoveSet from binaryCache" << endl;
    gettimeofday(&start, NULL);
    set1 = new NuPairMoveSetLibrary(libType, true,1);
    gettimeofday(&end, NULL);
    timeBinaryCache = end.tv_sec - start.tv_sec + (double)(end.tv_usec-start.tv_usec) /1e6;
    cout << "time consumed: " << timeBinaryCache << " seconds" << endl;
    cout << "first int in each vector of nbMoveList[6][9] moveIndexList:" << endl;
    for(int i=0; i<20; i++) {
        cout << set1->getMoveSet(6,9,1)->moveIndexList[i][0] << " ";
    }
    cout << endl;
    cout << "first int in each vector of revNbMoveList[6][9] moveIndexList:" << endl;
    for(int i=0; i<20; i++) {
        cout << set1->getMoveSet(6,9,-1)->moveIndexList[i][0] << " ";
    }
    cout << endl;
    cout << "first int in each vector of nnbMoveList[6][9] moveIndexList:" << endl;
    for(int i=0; i<20; i++) {
        cout << set1->getMoveSet(6,9,2)->moveIndexList[i][0] << " ";
    }
    cout << endl;
    cout << "oi sphereKeyMap1000 first 5 key-val pairs: " << endl;
    ii = 0;
	for(auto& it:set1->oi->getSKM1000()){
		cout << it.first << ":" << it.second << endl;
		ii++;
        if(ii>4) break;
	}
    set1->printMoveLibInfo();
    cout << "Reading BinaryTable is " << timeTxt/timeBinaryCache << "x boosting than reading txt." << endl;
    
    delete set1;
}