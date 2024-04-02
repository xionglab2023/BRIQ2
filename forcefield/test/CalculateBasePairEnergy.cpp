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



    char xx[200];
    string augc = "AUGC";

    cout << "load atLib" << endl;
    AtomLib* atLib = new AtomLib();

    cout << "load et1" << endl;
    ForceFieldPara* para1 = new ForceFieldPara();
    para1->bwTag = "adj";

    BasePair6DEnergyTable* et1 = new BasePair6DEnergyTable(para1,true, 1);
    et1->load(para1);


    for(int typeA=0;typeA<4;typeA++){
        for(int typeB=0;typeB<4;typeB++){
            sprintf(xx, "%c%c", augc[typeA], augc[typeB]);
            string cmListFile = "/public/home/pengx/briqx/basePair/nbPair/cluster/" + string(xx)+ ".cluster.cm";
            string output = "/public/home/pengx/briqx/basePair/nbPair/cluster/" + string(xx)+ ".cluster.ene";
                vector<CsMove> cmList;
    ifstream file;
    file.open(cmListFile.c_str(), ios::in);
    string line;
    while(getline(file, line)){
        CsMove cm(line);
        cmList.push_back(cm);
    }

    
    BaseRotamer* rotA = new BaseRotamer(typeA, atLib);
    BaseRotamer* rotB = new BaseRotamer(typeB, atLib);

    BasePairLib* bpLib = new BasePairLib();

   
    ofstream out;

    out.open(output.c_str(), ios::out);

    cout << "calculate energy: " << endl;

    for(int i=0;i<cmList.size();i++){
            CsMove cm = cmList[i];
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

        
            double e1 = et1->getEnergy(csA, csB, typeA, typeB, 2, minD);

            
            sprintf(xx, "%-8.4f", e1);
            out << string(xx) << endl;
        }
        out.close();
        delete rotA;
        delete rotB;
        }  
    }

    delete et1;
    delete para1;
    delete atLib;

}