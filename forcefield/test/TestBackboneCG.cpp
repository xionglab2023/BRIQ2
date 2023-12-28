
#include "geometry/localframe.h"
#include "geometry/TransMatrix.h"
#include "model/BaseDistanceMatrix.h"
#include "model/BasePairComposition.h"
#include "model/RiboseRotamerLib.h"
#include "forcefield/XPara.h"
#include "model/AtomLib.h"
#include "predNA/BRNode.h"
#include <time.h>
#include <stdio.h>

#include "model/StructureModel.h"
#include "forcefield/AtomicClashEnergy.h"
#include "forcefield/HbondEnergy.h"
#include "forcefield/BackboneConnectionEnergyCG.h"

using namespace NSPmodel;
using namespace NSPforcefield;


int main(int argc, char** argv){
	clock_t start = clock();
    ForceFieldPara* ffp = new ForceFieldPara();
    BackboneConnectionEnergyCG* bbcg = new BackboneConnectionEnergyCG(ffp);

    double ang1 = 120.0;
    double ang2 = 120.1;
    double dihed = 5.0;

    for(double d=2.0;d<6.0;d+=0.1){
        double e = bbcg->getEnergy(d, ang1, ang2, dihed);
        cout << d << " " << e << endl;
    }
    
}