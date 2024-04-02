/*
 * EnergyTable.cpp
 *
 */

#include <forcefield/RnaEnergyTable.h>

namespace NSPforcefield {

	RnaEnergyTable::~RnaEnergyTable(){
		if(this->acET != NULL) delete this->acET;
		if(this->bpET != NULL) delete this->bpET;
		if(this->pb != NULL) delete this->pb;
		if(this->roET != NULL) delete this->roET;
		if(this->hbET != NULL) delete this->hbET;
		if(this->bpcgET != NULL) delete this->bpcgET;
		if(this->bbcgET != NULL) delete this->bbcgET;
		if(this->bpLib != NULL) delete this->bpLib;
		
		if(deleteTag)
			delete this->para;

	}

} /* namespace NSPtest */
