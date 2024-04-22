
#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include "model/StructureModel.h"
#include "predNA/NuGraph.h"
#include "tools/CmdArgs.h"
#include "tools/ThreadPool.h"
#include "predNA/NuSampling.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredNA;
using namespace std;
using namespace NSPtools;
using namespace NSPthread;


int main(int argc, char** argv){

    BasePairLib* bpLib = new BasePairLib();
    EdgeInformationLib* eiLib = new EdgeInformationLib(bpLib);
    vector<string> list = eiLib->keyList;
    for(int i=0;i<list.size();i++){
        string key = list[i];
        int num = eiLib->eiMap[key]->validClusterNum;
        double weight = eiLib->eiMap[key]->weight;
        cout << key  << " " << num << " " << weight << endl;

        
    }
}