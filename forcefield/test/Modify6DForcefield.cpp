#include <vector>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include "geometry/OrientationIndex.h"
#include "model/PolarGroupLib.h"
#include "forcefield/HbondEnergy.h"
#include "forcefield/AtomicClashEnergy.h"
#include "forcefield/Xtb6dEnergy.h"
#include "forcefield/BasePair6DEnergyTable.h"
#include "model/BaseRotamer.h"
#include "model/BasePairLib.h"

using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;
using namespace NSPforcefield;

int main(int argc, char** argv){

    string listFile = string(argv[1]);
    ifstream in;
    in.open(listFile.c_str(), ios::in);
    string line;
    string pairType;
    getline(in, line);
    pairType = line;
    set<int> selectClusters;
    while(getline(in, line)){
        selectClusters.insert(atoi(line.c_str()));
    }
    in.close();

    cout << pairType << endl;
    map<char, int> augc;
    augc['A'] = 0;
    augc['U'] = 1;
    augc['G'] = 2;
    augc['C'] = 3;
    int typeA = augc[pairType[0]];
    int typeB = augc[pairType[1]];

    string outFile = "/public/home/pengx/cpp/briqx/data/pairEne/nnb/clustersNeedToOpt/"+pairType+".ene-adj";
    ofstream out;
    out.open(outFile.c_str(), ios::out);

    set<int>::iterator it;
    for(it = selectClusters.begin();it!=selectClusters.end();++it){
        cout << *it << endl;
    }

    ForceFieldPara* para = new ForceFieldPara();
    para->bwTag = "raw";
    BasePair6DEnergyTable* et = new BasePair6DEnergyTable(para, true, 1);
	et->load(para);

    string eneFile = "/public/home/pengx/cpp/briqx/data/pairEne/nnb/oldAdj/"+pairType+".ene-adj";
    in.open(eneFile, ios::in);
    int idA, idB, clusterID;
    double ene, ene2;
    OrientationIndex oi;

    while(in >> idA >> idB >> ene >> clusterID){
        if(selectClusters.find(clusterID) != selectClusters.end()){
            CsMove cm = oi.index2000ToCsMove(idA, idB);
            LocalFrame cs1;
            LocalFrame cs2 = cs1 + cm;
            ene2 = et->getEnergy(cs1, cs2, typeA, typeB, 2, 4.0) / para->wtBp2 * 1.34;
            if(ene2 < -0.01)
                out << idA << " " << idB  << " " << ene2 << " " << clusterID << endl;
        }
        else {
             out << idA << " " << idB  << " " << ene << " " << clusterID << endl;
        }
    }
    in.close();
    out.close();

}