/*
 * RnsEnergyTableSimple.h
 *
 */

#ifndef FORCEFIELD_RNAENERGYTABLESIMPLE_H_
#define FORCEFIELD_RNAENERGYTABLESIMPLE_H_

#include <time.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include "model/BaseDistanceMatrix.h"
#include "dataio/datapaths.h"
#include "forcefield/XPara.h"
#include "forcefield/AtomicClashEnergy.h"
#include "forcefield/BasePair6DEnergyTable.h"
#include "forcefield/RiboseOxygenEnergyTable.h"
#include "forcefield/ForceFieldPara.h"
#include "RnaAtomicEnergyTable.h"

namespace NSPforcefield {

using namespace std;
using namespace NSPmodel;


class RnaEnergyTableSimple{
public:

	AtomicClashEnergy* acET;
	RiboseOxygenEnergyTable* roET;
	XPara para;


	RnaEnergyTableSimple(){
		ForceFieldPara* ffp = new ForceFieldPara();
		acET = new AtomicClashEnergy(ffp);
		roET = new RiboseOxygenEnergyTable();
	}


	virtual ~RnaEnergyTableSimple();
};

} /* namespace NSPforceField */

#endif /* FORCEFIELD_RNAENERGYTABLESIMPLE_H_ */
