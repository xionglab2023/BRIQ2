#include <vector>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include "geometry/OrientationIndex.h"
#include "model/PolarGroupLib.h"
#include "forcefield/HbondEnergy.h"
#include "forcefield/AtomicClashEnergy.h"
#include "forcefield/Xtb6dEnergy.h"
#include "model/BaseRotamer.h"

using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;
using namespace NSPforcefield;


int main(int argc, char** argv){

    /*
     *  usage: xtbEnergy $GROUP_TYPE_A $GROUP_TYPE_B $OI_TYPE $OI_INDEX ($OI_INDEX2) 
     */

    int typeA = atoi(argv[1]);
    int typeB = atoi(argv[2]);
    int idTag = atoi(argv[3]);


    string outfile = string(argv[4]);

    char xx[100];
    OrientationIndex x;

    ofstream out;
    out.open(outfile, ios::out);
    int i,j,k,l;
    CsMove cm;
    LocalFrame csA;
    LocalFrame csB;
    AtomLib* atLib = new AtomLib();
    BaseRotamer* rotA = new BaseRotamer(typeA, atLib);
    BaseRotamer* rotB = new BaseRotamer(typeB, atLib);
    BaseConformer* confA;
    BaseConformer* confB;
    ForceFieldPara* para = new ForceFieldPara();
    Xtb6dEnergy* et = new Xtb6dEnergy(para);

    OrientationIndex oi;
    double minD = 99.9;
    double d;
    double ene;
    int index1, index2;
    for(index1=0;index1<2250;index1++) {
        cout << index1 << endl;
        if(index1 % 50 != idTag) continue;

        for(index2=0;index2<4000000;index2++){
            cm = oi.index2000ToCsMove(index1, index2);
            int id500 = oi.moveToIndex500(cm);
            csB = csA + cm;
            confA = new BaseConformer(rotA, csA);
            confB = new BaseConformer(rotB, csB);
            minD = 99.9;
            for(i=0;i<confA->rot->atomNum;i++){
                for(j=0;j<confB->rot->atomNum;j++){
                    d = confA->coords[i].distance(confB->coords[j]);
                    if(d < minD){
                     minD = d;
                    }
                }
            }

            if(typeA <= typeB)
                ene = et->getEnergy(csA, csB, typeA, typeB, minD);
            else
                ene = et->getEnergy(csB, csA, typeB, typeA, minD);

            if(ene < -2.5) {
                out << index1 << " " << index2 << " " << ene << endl;
            }
            delete confA;
            delete confB;
        }
    }

    out.close();

    delete atLib;
    delete rotA;
    delete rotB;
    delete para;
    delete et;
}