/*
 *  BackboneConnectionEnergyCG.h
 * Created on 2023/12/26
 * Author: pengx
 */

#ifndef FORCEFIELD_BACKBONECONNECTIONENERGYCG_H_
#define FORCEFIELD_BACKBONECONNECTIONENERGYCG_H_

#include "dataio/datapaths.h"
#include "forcefield/ForceFieldPara.h"
#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include "forcefield/BasePair6DEnergyTable.h"
#include <time.h>
#include <math.h>
#include <map>
#include <fstream>
#include <sstream>

namespace NSPforcefield {

using namespace std;
using namespace NSPgeometry;

class BackboneConnectionEnergyCG{
public:

    /* 
     * 15*40*40*120 = 2880000
     * C1'A - O3'A - C5'B - C1'B
     * d: O3' - C5', range 2.6~4.1
     * ang1: C1'-O3'-C5', range 60~180
     * ang2: O3'-C5'-C1', range 60~180
     * dihed: C1'-O3'-C5'-C1', range -180~180
     */

    double et[2880000];
    BackboneConnectionEnergyCG(ForceFieldPara* pa);
    double getEnergy(double d, double ang1, double ang2, double dihed);
    virtual ~BackboneConnectionEnergyCG();
};

}

#endif