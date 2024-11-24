
#include <predNA/EdgeClusterRegister.h>

namespace NSPpredNA {

EdgeClusterRegister::EdgeClusterRegister(int baseTypeA, int baseTypeB, BasePairLib* bpLib){


    this->typeA = baseTypeA%4;
    this->typeB = baseTypeB%4;

    int pairType = typeA*4+typeB;
    this->totalClusterNum = bpLib->nnbBasePairNum[typeA*4+typeB];

    for(int i=0;i<totalClusterNum;i++){
        if(bpLib->nnbEnergyWithOxy[pairType][i] < -2.0) {
            this->lowEnergyClusterIDList.push_back(i);
        }
    }

    this->clusterIDCount = new int[totalClusterNum];
    this->pCluster = new double[totalClusterNum];

    for(int i=0;i<totalClusterNum;i++){
        this->clusterIDCount[i] = 0;
        this->pCluster[i] = 0.0;
    }
    this->totalCount = 0;
}

void EdgeClusterRegister::clear(){
    for(int i=0;i<totalClusterNum;i++){
        this->clusterIDCount[i] = 0;
        this->pCluster[i] = 0.0;
    }
    this->totalCount = 0;
}

void EdgeClusterRegister::record(int clusterID){

    clusterIDCount[clusterID]++;
    totalCount ++;
}


void EdgeClusterRegister::print(){

    int n1, n2;

    if(totalCount == 0)
        return;

    for(int i=0;i<totalClusterNum;i++){
        n1 = this->clusterIDCount[i];
        if(n1 > 0 || n2 > 0)
            cout << "clusterID: " << i << " " << n1 << endl;
    }
}

EdgeClusterRegister::~EdgeClusterRegister(){
    delete [] clusterIDCount;
    delete [] pCluster;
}

}