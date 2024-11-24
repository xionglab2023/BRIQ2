#ifndef PREDNA_EDGECLUSTERREGISTER_H_
#define PREDNA_EDGECLUSTERREGISTER_H_

#include "model/BasePairLib.h"


namespace NSPpredNA {

using namespace NSPmodel;

class EdgeClusterRegister{
public:

    int typeA;
    int typeB;
    int totalClusterNum;
    vector<int> lowEnergyClusterIDList;

    int* clusterIDCount;

    int totalCount;

    double* pCluster;



    EdgeClusterRegister(int baseTypeA, int baseTypeB, BasePairLib* bpLib);
    
    void clear();
   
    void record(int clusterID);
    void print();
    virtual ~EdgeClusterRegister();
};

} /* namespace NSPpredNA */

#endif /* PREDNA_EDGEINFORMATION_H_ */