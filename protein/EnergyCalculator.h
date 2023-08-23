/*
 * EnergyCalculator.h
 *
 */

#ifndef PREDPROT_ENERGYCALCULATOR_H_
#define PREDPROT_ENERGYCALCULATOR_H_

#include "model/ResConformer.h"
#include "forcefield/ProEnergyTable.h"
#include "para/ProParameter.h"
#include "para/DesignPara.h"

namespace NSPprotein {

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpara;
using namespace std;

class EnergyCalculator {

public:
	ProEnergyTable* et;
	DesignPara* dp;
	ProParameter* pa;


	EnergyCalculator(ProEnergyTable* et, DesignPara* dp, ProParameter* pa){
		this->et = et;
		this->pa = pa;
		this->dp = dp;
	}

	double getEnergyDesign(ResConformer* confA, ResConformer* confB, float meanSai, int sep){
		return 0.0;
	}

	double getEnergyABACUS(ResConformer* confA, ResConformer* confB, float meanSai, int sep){
		return 0.0;
	}

	double getEnergyBbScDesign(ResConformer* confA, ResConformer* confB, float meanSai, int sep){
		return 0.0;
	}

	double getEnergyScSc(ResConformer* scConfA, ResConformer* scConfB, int sep){
		return 0.0;
	}
	double getEnergyBbSc(ResConformer* bbConfA, ResConformer* scConfB, int sep){
		return 0.0;
	}

//	double getEnergyScScDesign(ResScConformer* scConfA, ResScConformer* scConfB, float meanSai, int sep);
//	double getEnergyBbScDesign(ResBBConformer* bbConfA, ResScConformer* scConfB, float meanSai, int sep);

	double printEnergyDesign(int idA, int idB, ResConformer* confA, ResConformer* confB, float meanSai, int sep){
		return 0.0;
	}

//	void printEnergyScScDesign(int idA, int idB, ResScConformer* scConfA, ResScConformer* scConfB, float meanSai, int sep, ProteinAtomLib* atLib);
//	void printEnergyBbScDesign(int idA, int idB, ResBBConformer* bbConfA, ResScConformer* scConfB, float meanSai, int sep, ProteinAtomLib* atLib);

	double getEnergyBbBb(ResBBConformer* bbConfA, ResBBConformer* bbConfB, int sep){
		return 0.0;
	}
	double getEnergyScScPolar(ResScConformer* scConfA, ResScConformer* scConfB, int sep){
		return 0.0;
	}
	double getEnergyBbScPolar(ResBBConformer* bbConfA, ResScConformer* scConfB, int sep){
		return 0.0;
	}

	virtual ~EnergyCalculator();

};


} /* namespace NSPforcefield */

#endif /* PREDPROT_ENERGYCALCULATOR_H_ */
