/*
 * EnergyArray.h
 *
 *  Created on: 2022Äê7ÔÂ16ÈÕ
 *      Author: pengx
 */

#ifndef PREDPROT_ENERGYARRAY_H_
#define PREDPROT_ENERGYARRAY_H_

#include <iostream>

namespace NSPprotein {
using namespace std;


class EnergyArray {

public:
	float* ea;
	int choiceNum;
	EnergyArray(int choiceNum);
	int getChoiceNum();
	void setEnergy(int choice, float e);
	float getEnergy(int choice);
	virtual ~EnergyArray();
};

} /* namespace NSPprotein */

#endif /* PREDPROT_ENERGYARRAY_H_ */
