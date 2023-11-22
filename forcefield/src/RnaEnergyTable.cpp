/*
 * EnergyTable.cpp
 *
 */

#include <forcefield/RnaEnergyTable.h>

namespace NSPforcefield {

	RnaEnergyTable::~RnaEnergyTable(){
		delete this->acET;
		delete this->bpET;
		delete this->pb;
		delete this->roET;
		delete this->hbET;
		if(deleteTag)
			delete this->para;

	}

} /* namespace NSPtest */
