/*
 * ProEnergyTable.h
 *
 */

#ifndef FORCEFIELD_PROENERGYTABLE_H_
#define FORCEFIELD_PROENERGYTABLE_H_

#include "forcefield/ProAtomicEnergyTable.h"
#include "para/ProParameter.h"
#include "para/DesignPara.h"

namespace NSPforcefield {

using namespace NSPpara;

class ProEnergyTable {
public:

	ProAtomicEnergyTable* etAt;

	ProEnergyTable();
	ProEnergyTable(DesignPara* dp);
	virtual ~ProEnergyTable();
};

} /* namespace NSPforcefield */

#endif /* FORCEFIELD_PROENERGYTABLE_H_ */
