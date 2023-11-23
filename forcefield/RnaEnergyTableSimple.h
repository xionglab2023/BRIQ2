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
#include "forcefield/PO3Builder.h"
#include "forcefield/AtomicClashEnergy.h"
#include "forcefield/BasePair6DEnergyTable.h"
#include "forcefield/RiboseOxygenEnergyTable.h"
#include "forcefield/ForceFieldPara.h"
#include "RnaAtomicEnergyTable.h"
#include "forcefield/HbondEnergy.h"

namespace NSPforcefield {

using namespace std;
using namespace NSPmodel;


class RnaEnergyTableSimple{
public:

	PO3Builder* pb;
	AtomicClashEnergy* acET;
	RiboseOxygenEnergyTable* roET;
	HbondEnergy* hbET;
	ForceFieldPara* para;


	RnaEnergyTableSimple(ForceFieldPara* para){

		this->para = para;
		pb = new PO3Builder(para);
		acET = new AtomicClashEnergy(para);
		roET = new RiboseOxygenEnergyTable(para);
		hbET = new HbondEnergy(para);
	}


	virtual ~RnaEnergyTableSimple();
};

} /* namespace NSPforceField */

#endif /* FORCEFIELD_RNAENERGYTABLESIMPLE_H_ */
