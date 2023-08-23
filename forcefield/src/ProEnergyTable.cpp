/*
 * ProEnergyTable.cpp
 *
 */

#include "forcefield/ProEnergyTable.h"

namespace NSPforcefield {

ProEnergyTable::ProEnergyTable() {
	this->etAt = new ProAtomicEnergyTable();

}

ProEnergyTable::ProEnergyTable(DesignPara* dp){
	this->etAt = new ProAtomicEnergyTable(dp);
}

ProEnergyTable::~ProEnergyTable() {
	// TODO Auto-generated destructor stub
	delete etAt;
}

} /* namespace NSPforcefield */
