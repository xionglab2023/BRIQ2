/*
 * RotamerLib.cpp
 *
 *  Created on: 2023Äê8ÔÂ14ÈÕ
 *      Author: nuc
 */

#include <model/RotamerLib.h>

namespace NSPmodel {


RotamerLib::~RotamerLib() {
	// TODO Auto-generated destructor stub
	if(this->bbRotLib != NULL)
		delete this->bbRotLib;
	if(this->scRotLib != NULL)
		delete this->scRotLib;
	if(this->pepRotLib != NULL)
		delete this->pepRotLib;
	if(this->baseRotLib != NULL)
		delete this->baseRotLib;
	if(this->riboseRotLib != NULL)
		delete this->riboseRotLib;
	if(this->phoRotLib != NULL)
		delete this->phoRotLib;
}

} /* namespace NSPmodel */
