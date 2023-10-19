/*
 * EnergyTable.h
 *
 */

#ifndef FORCEFIELD_RNAENERGYTABLE_H_
#define FORCEFIELD_RNAENERGYTABLE_H_

#include "model/BaseDistanceMatrix.h"
#include "dataio/datapaths.h"
#include "forcefield/XPara.h"

#include "forcefield/BasePair6DEnergyTable.h"
#include "forcefield/PO3Builder.h"
#include "forcefield/RiboseOxygenEnergyTable.h"
#include "forcefield/AtomicClashEnergy.h"
#include "forcefield/HbondEnergy.h"

#include <time.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include "RnaAtomicEnergyTable.h"

namespace NSPforcefield {

using namespace std;
using namespace NSPmodel;


class RnaEnergyTable{
public:
	BasePair6DEnergyTable* bpET;
	AtomicClashEnergy* acET;
	PO3Builder* pb;
	RiboseOxygenEnergyTable* roET;
	HbondEnergy* hbET;
	ForceFieldPara* para;

	RnaEnergyTable(){
		this->para = new ForceFieldPara();
		bpET = new BasePair6DEnergyTable(para);
		acET = new AtomicClashEnergy(para);
		roET = new RiboseOxygenEnergyTable(para);
		pb = new PO3Builder(para);
		hbET = new HbondEnergy(para);
	}

	RnaEnergyTable(const string& paraFile){

		this->para = new ForceFieldPara();
		bpET = new BasePair6DEnergyTable(para);
		acET = new AtomicClashEnergy(para);
		roET = new RiboseOxygenEnergyTable(para);
		pb = new PO3Builder(para);
		hbET = new HbondEnergy(para);
	}

	virtual ~RnaEnergyTable();
};

} /* namespace NSPforceField */

#endif /* FORCEFIELD_RNAENERGYTABLE_H_ */
