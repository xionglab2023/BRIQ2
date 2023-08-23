/*
 * EnergyTable2.cpp
 *
 */

#include <forcefield/RnaEnergyTableSimple.h>

namespace NSPforcefield {

RnaEnergyTableSimple::~RnaEnergyTableSimple() {
	// TODO Auto-generated destructor stub

	delete this->acET;
	delete this->roET;
}

} /* namespace NSPforcefield */
