/*
 * RotamerLib.h
 *
 *  Created on: 2023Äê8ÔÂ14ÈÕ
 *      Author: nuc
 */

#ifndef MODEL_ROTAMERLIB_H_
#define MODEL_ROTAMERLIB_H_

#include "model/ResBBRotamerLib.h"
#include "model/ResScRotamerLib.h"
#include "model/ResPeptideRotamerLib.h"
#include "model/BaseRotamerLib.h"
#include "model/RiboseRotamerLib.h"
#include "model/PhosphateRotamerLib.h"
#include "forcefield/ForceFieldPara.h"

namespace NSPmodel {

class RotamerLib {
public:

	ResBBRotamerLib* bbRotLib;
	ResScRotamerLib* scRotLib;
	ResPeptideRotamerLib* pepRotLib;
	BaseRotamerLib* baseRotLib;
	RiboseRotamerLib* riboseRotLib;
	PhosphateRotamerLib* phoRotLib;

	RotamerLib(){
		AtomLib* atLib = new AtomLib();
		this->baseRotLib = new BaseRotamerLib(atLib);
		this->riboseRotLib = new RiboseRotamerLib();
		this->phoRotLib = new PhosphateRotamerLib();

		this->bbRotLib = NULL;
		this->scRotLib = NULL;
		this->pepRotLib = NULL;
	}

	RotamerLib(ForceFieldPara* para){
		AtomLib* atLib = new AtomLib();
		this->baseRotLib = new BaseRotamerLib(atLib);
		this->riboseRotLib = new RiboseRotamerLib(para);
		this->phoRotLib = new PhosphateRotamerLib();

		this->bbRotLib = NULL;
		this->scRotLib = NULL;
		this->pepRotLib = NULL;
	}

	RotamerLib(const string& tag){
		if(tag == "rna"){
			AtomLib* atLib = new AtomLib();
			this->baseRotLib = new BaseRotamerLib(atLib);
			this->riboseRotLib = new RiboseRotamerLib();
			this->phoRotLib = new PhosphateRotamerLib();
			this->bbRotLib = NULL;
			this->scRotLib = NULL;
			this->pepRotLib = NULL;
		}
		else if(tag == "protein"){
			this->bbRotLib = new ResBBRotamerLib();
			this->scRotLib = new ResScRotamerLib();
			this->pepRotLib = new ResPeptideRotamerLib();
			this->baseRotLib = NULL;
			this->riboseRotLib = NULL;
			this->phoRotLib = NULL;
		}
		else if(tag == "all") {
			AtomLib* atLib = new AtomLib();
			this->bbRotLib = new ResBBRotamerLib();
			this->scRotLib = new ResScRotamerLib();
			this->pepRotLib = new ResPeptideRotamerLib();
			this->baseRotLib = new BaseRotamerLib(atLib);
			this->riboseRotLib = new RiboseRotamerLib();
			this->phoRotLib = new PhosphateRotamerLib();
		}
		else {
			AtomLib* atLib = new AtomLib();
			this->bbRotLib = new ResBBRotamerLib();
			this->scRotLib = new ResScRotamerLib();
			this->pepRotLib = new ResPeptideRotamerLib();
			this->baseRotLib = new BaseRotamerLib(atLib);
			this->riboseRotLib = new RiboseRotamerLib();
			this->phoRotLib = new PhosphateRotamerLib();
		}
	}

	virtual ~RotamerLib();
};

} /* namespace NSPmodel */

#endif /* MODEL_ROTAMERLIB_H_ */
