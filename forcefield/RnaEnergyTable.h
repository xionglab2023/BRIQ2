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
#include "forcefield/BasePair6DEnergyTableCG.h"
#include "forcefield/PO3Builder.h"
#include "forcefield/RiboseOxygenEnergyTable.h"
#include "forcefield/AtomicClashEnergy.h"
#include "forcefield/HbondEnergy.h"
#include "forcefield/BackboneConnectionEnergyCG.h"

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
	BasePair6DEnergyTableCG* bpcgET;
	AtomicClashEnergy* acET;
	PO3Builder* pb;
	RiboseOxygenEnergyTable* roET;
	HbondEnergy* hbET;
	BackboneConnectionEnergyCG* bbcgET;

	ForceFieldPara* para;
	bool deleteTag;

	RnaEnergyTable(){
		this->para = new ForceFieldPara();
		deleteTag = true;
		bpET = new BasePair6DEnergyTable(para);
		bpcgET = new BasePair6DEnergyTableCG(para);
		acET = new AtomicClashEnergy(para);
		roET = new RiboseOxygenEnergyTable(para);
		pb = new PO3Builder(para);
		hbET = new HbondEnergy(para);
		bbcgET = new BackboneConnectionEnergyCG(para);
	}

	RnaEnergyTable(const string& paraFile){

		this->para = new ForceFieldPara();
		deleteTag = true;
		bpET = new BasePair6DEnergyTable(para);
		bpcgET = new BasePair6DEnergyTableCG(para);
		acET = new AtomicClashEnergy(para);
		roET = new RiboseOxygenEnergyTable(para);
		pb = new PO3Builder(para);
		hbET = new HbondEnergy(para);
		bbcgET = new BackboneConnectionEnergyCG(para);
	}

	RnaEnergyTable(ForceFieldPara* para){
		this->para = para;
		deleteTag = false;
		bpET = new BasePair6DEnergyTable(para);
		bpcgET = new BasePair6DEnergyTableCG(para);
		acET = new AtomicClashEnergy(para);
		roET = new RiboseOxygenEnergyTable(para);
		pb = new PO3Builder(para);
		hbET = new HbondEnergy(para);
		bbcgET = new BackboneConnectionEnergyCG(para);
	}

	virtual ~RnaEnergyTable();
};

} /* namespace NSPforceField */

#endif /* FORCEFIELD_RNAENERGYTABLE_H_ */
