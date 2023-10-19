/*
 * EnergyTable2.cpp
 *
 */

#include <forcefield/RnaEnergyTableSimple.h>

namespace NSPforcefield {

RnaEnergyTableSimple::~RnaEnergyTableSimple() {
	// TODO Auto-generated destructor stub

	delete pb;
	delete hbET;
	delete acET;
	delete roET;
}

} /* namespace NSPforcefield */
