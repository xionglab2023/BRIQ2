/*
 * EnergyTable.h
 *
 */

#ifndef FORCEFIELD_RNAENERGYTABLE_H_
#define FORCEFIELD_RNAENERGYTABLE_H_

#include "model/BaseDistanceMatrix.h"
#include "model/BasePairLib.h"
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

		bpET = NULL;
		bpcgET = NULL;
		acET = NULL;
		roET = NULL;
		pb = NULL;
		hbET = NULL;
		bbcgET = NULL;
	}

	RnaEnergyTable(const string& paraFile){

		this->para = new ForceFieldPara();
		deleteTag = true;
		bpET = NULL;
		bpcgET = NULL;
		acET = NULL;
		roET = NULL;
		pb = NULL;
		hbET = NULL;
		bbcgET = NULL;
	}

	RnaEnergyTable(ForceFieldPara* para){
		this->para = para;
		deleteTag = false;
		bpET = NULL;
		bpcgET = NULL;
		acET = NULL;
		roET = NULL;
		pb = NULL;
		hbET = NULL;
		bbcgET = NULL;
	}

	void loadEnergyWithout6D(){

		acET = new AtomicClashEnergy(para);
		roET = new RiboseOxygenEnergyTable(para);
		pb = new PO3Builder(para);
		hbET = new HbondEnergy(para);	
	}

	void loadAtomicEnergy(){
		cout << "load bp6d" << endl;
		bpET = new BasePair6DEnergyTable(para, true, 1);
		bpET->load(para);
		cout << "load clash" << endl;
		acET = new AtomicClashEnergy(para);
		cout << "load ribose oxy" << endl;
		roET = new RiboseOxygenEnergyTable(para);
		cout << "load pb" << endl;
		pb = new PO3Builder(para);
		cout << "load hb" << endl;
		hbET = new HbondEnergy(para);
	}

	void updateAtomic(ForceFieldPara* para){
		if(this->acET != NULL)
			delete this->acET;
		this->acET = new AtomicClashEnergy(para);

		if(this->bpET != NULL)
			delete hbET;
		this->hbET = new HbondEnergy(para);

	}

	void loadCoarseGrainedEnergy(){
		pb = new PO3Builder(para);
		bpcgET = new BasePair6DEnergyTableCG(para, true, 1);
		bpcgET->load();
		acET = new AtomicClashEnergy(para);
		bbcgET = new BackboneConnectionEnergyCG(para);
	}	

	virtual ~RnaEnergyTable();
};

} /* namespace NSPforceField */

#endif /* FORCEFIELD_RNAENERGYTABLE_H_ */
