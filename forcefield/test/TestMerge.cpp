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

    string outfile = string(argv[3]);

    MergeThreeNnbEnergyTable* et = new MergeThreeNnbEnergyTable(typeA, typeB);
    et->mergeEnergy(outfile);
    delete et;
    
}