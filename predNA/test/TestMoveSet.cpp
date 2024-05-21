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


    struct timeval start, end;
    double timeTxt, timeBinaryTable, timeBinaryCache;

    string libType = string(argv[1]);

	NuPairMoveSetLibrary* moveLib = new NuPairMoveSetLibrary(libType, true, 1);
	moveLib->load();

    BasePairLib* bpLib = new BasePairLib(libType);

    string augc= "AUGC";
    OrientationIndex* oi = new OrientationIndex();
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            string pairType = augc.substr(i, 1) + augc.substr(j, 1);
            int n = bpLib->nnbBasePairNum[i*4+j];
            for(int k=0;k<n;k++){
                double minD = 99.9;
                double maxD = 0.0;
                int moveNum = 0;
                BaseDistanceMatrix dm1 = bpLib->nnbDMClusterCenters[i*4+j][k];
                IndividualNuPairMoveSet* set = moveLib->getMoveSet(i*4+j, k, 2);
                for(int x=0;x<50;x++){
                    int n2 = set->moveIndexList[x].size();
                    if(n2 == 0) {
                        printf("%s %4d\n", pairType.c_str(), k); 
                    }
                    moveNum += n2;
                    for(int y=0;y<n2;y++){
                        int index = set->moveIndexList[x][y];
                        CsMove cm = oi->index1000ToCsMove(index);
                        BaseDistanceMatrix dm2(cm);
                        double d = dm1.distanceTo(dm2);
                        if(d < minD){
                            minD = d;
                        }
                        if(d > maxD){
                            maxD = d;
                        }
                    }
                }
                printf("%s %4d %6.3f %6.3f %7d\n", pairType.c_str(), k, minD, maxD, moveNum);
            }
        }
    }

}