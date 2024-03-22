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

    int typeA = atoi(argv[1]);
    int typeB = atoi(argv[2]);
    string nbnnb = string(argv[3]);

    char xx[200];
    string augc = "AUGC";

    cout << "load atLib" << endl;
    AtomLib* atLib = new AtomLib();

    cout << "load et1" << endl;
    ForceFieldPara* para1 = new ForceFieldPara();
    para1->bwTag = "raw";
    BasePair6DEnergyTable* et1 = new BasePair6DEnergyTable(para1,true, 1);
    et1->load(para1);

    cout << "load et2" << endl;
    ForceFieldPara* para2 = new ForceFieldPara();
    para2->bwTag = "adj";
    BasePair6DEnergyTable* et2 = new BasePair6DEnergyTable(para2, true, 1);
    et2->load(para2);

    
    OrientationIndex* oi = new OrientationIndex();
    
    BaseRotamer* rotA = new BaseRotamer(typeA, atLib);
    BaseRotamer* rotB = new BaseRotamer(typeB, atLib);

    BasePairLib* bpLib = new BasePairLib();

    int clusterNum = bpLib->nnbBasePairNum[typeA*4+typeB];
    if(nbnnb == "nb") {
        clusterNum = bpLib->nbBasePairNum[typeA*4+typeB];
    }

    delete bpLib;
    ifstream file;
    ofstream out;

    cout << "calculate energy: " << endl;

    for(int clusterID=0;clusterID<clusterNum;clusterID++){
        vector<CsMove> moveList;
        sprintf(xx, "/public/home/pengx/cpp/briqx/data/pairMove4/%s/%c%c%d.move", nbnnb.c_str(), augc[typeA], augc[typeB], clusterID);
        string moveFile = string(xx);

        sprintf(xx, "/public/home/pengx/briqx/bpEnergyAnalysis/%s/%c%c/%c%c%d.ene", nbnnb.c_str(), augc[typeA], augc[typeB],augc[typeA], augc[typeB], clusterID);
        string outFile = string(xx);

        out.open(outFile.c_str(), ios::out);

        file.open(moveFile.c_str(), ios::in);
        int id1000, idc;
        double p;
        while(file >> id1000 >> p >> idc){
            moveList.push_back(oi->index1000ToCsMove(id1000));
        }
        file.close();

        vector<double> e1List;
        vector<double> e2List;
        for(int i=0;i<moveList.size();i++){
            CsMove cm = moveList[i];
            LocalFrame csA;
            LocalFrame csB = csA.add(cm);

            BaseConformer* confA = new BaseConformer(rotA, csA);
            BaseConformer* confB = new BaseConformer(rotB, csB);
            double minD = 99.9;
            for(int j=0;j<rotA->atomNum;j++){
                for(int k=0;k<rotB->atomNum;k++){
                    double d = confA->coords[j].distance(confB->coords[k]);
                    if(d < minD){
                        minD = d;
                    }
                }   
            }
            delete confA;
            delete confB;

            int sep = 2;
            if(nbnnb == "nb")
                sep = 1;

            double e1, e2;
            e1 = et1->getEnergy(csA, csB, typeA, typeB, sep, minD);

            if(sep == 1)
                e2 = et2->getEnergy(csA, csB, typeA, typeB, 2, minD);
            sprintf(xx, "%-8.4f %-8.4f", e1, e2);
            out << string(xx) << endl;
        }
        out.close();
    }


    delete rotA;
    delete rotB;
    delete oi;
    delete et1;
    delete et2;
    delete para1;
    delete para2;
    delete atLib;

}