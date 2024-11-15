
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

    this->clusterExistence = new int[totalClusterNum];
    this->clusterExistenceReverse = new int[totalClusterNum];

    for(int i=0;i<totalClusterNum;i++){
        this->clusterExistence[i] = 0;
    }
}

void EdgeClusterRegister::clear(){
    for(int i=0;i<totalClusterNum;i++){
        this->clusterExistence[i] = 0;
        this->clusterExistenceReverse[i] = 0;
    }
}

void EdgeClusterRegister::record(int clusterID){
    clusterExistence[clusterID]++;
}

void EdgeClusterRegister::recordReverse(int clusterID){
    clusterExistenceReverse[clusterID]++;
}

EdgeClusterRegister::~EdgeClusterRegister(){
    delete [] clusterExistence;
    delete [] clusterExistenceReverse;
}

}