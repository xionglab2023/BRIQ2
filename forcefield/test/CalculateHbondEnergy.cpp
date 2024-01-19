/*
 * CalculateHbondEnergy.cpp
 *
 *  Created on: 2024,1,9
 *      Author: pengx
 */

#include <vector>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include "geometry/OrientationIndex.h"
#include "model/PolarGroupLib.h"
#include "forcefield/HbondEnergy.h"
#include "forcefield/AtomicClashEnergy.h"
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
    int index1 = atoi(argv[3]);
    int index2;
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

    int npA = rotA->polarAtomNum;
    int npB = rotB->polarAtomNum;
    ForceFieldPara* para = new ForceFieldPara();
    HbondEnergy* hbet = new HbondEnergy(para);
    AtomicClashEnergy* acet = new AtomicClashEnergy(para);
    double dd;
    double e, totalE;
    int hbNum;
    double clashEnergy;

    for(i=0;i<2000;i++){
        for(j=0;j<2000;j++){
            index2 = i*2000+j;
            cm = x.index2000ToCsMove(index1, index2);
            csB = csA.add(cm);

            confA = new BaseConformer(rotA, csA);
            confB = new BaseConformer(rotB, csB);

            totalE = 0.0;
            clashEnergy = 0.0;
            hbNum = 0;
            for(k=0;k<npA;k++){
                for(l=0;l<npB;l++){
                    dd = confA->csPolar[k].origin_.squaredDistance(confB->csPolar[l].origin_);
                    if(dd > 13.0) continue;
                    e = hbet->getEnergy(confA->rot->polarAtomUniqueID[k], confA->csPolar[k], confB->rot->polarAtomUniqueID[l], confB->csPolar[l]);
                    totalE += e;
                    if(e < -0.1) {
                        hbNum++;
                    }
                }
            }

            for(k=0;k<confA->rot->atomNum;k++){
                for(l=0;l<confB->rot->atomNum;l++){
                    dd = confA->coords[k].squaredDistance(confB->coords[l]);
                    clashEnergy += acet->getBaseBaseEnergy(typeA, k, typeB, l, dd, 2);
                }
            }

            for(k=0;k<confA->rot->atomNum;k++){
                dd = confA->coords[k].squaredDistance(confB->cs1.origin_);
                clashEnergy += acet->getBaseRiboseEnergy(typeA, k, 0, dd, 2);
            }

            for(k=0;k<confB->rot->atomNum;k++){
                dd = confB->coords[k].squaredDistance(confA->cs1.origin_);
                clashEnergy += acet->getBaseRiboseEnergy(typeB, k, 0, dd, 2);
            }

            delete confA;
            delete confB;

            totalE = totalE+0.1;

            if(clashEnergy + totalE > 0.0) continue;

            if(totalE < 0){
                sprintf(xx, "%d %d %8.4f %d", index1, index2, totalE, hbNum);
                out << string(xx) << endl;
            }
        }
    }

    out.close();

    delete atLib;
    delete rotA;
    delete rotB;
    delete hbet;
    delete acet;
    delete para;
}


