/*
 * AAProbabilityArray.h
 *
 */

#ifndef MATH_AAPROBABILITYARRAY_H_
#define MATH_AAPROBABILITYARRAY_H_

#include "math/SubstitutionMatrix.h"
#include <math.h>
#include "model/ResName.h"

namespace NSPmath {

class AAProbabilityArray {
public:
	double pa[20];
	double sampelNum = -1;
	AAProbabilityArray();
	AAProbabilityArray(double* p, double N);
	AAProbabilityArray(int* count);
	void normalize();
	void addPseudoCount(SubstitutionMatrix* subM, int pseudoNum);
	double entropy();
	double distance(AAProbabilityArray* other);
	void printAAName();
	void print();

	virtual ~AAProbabilityArray();
};

} /* namespace NSPmath */

#endif /* MATH_AAPROBABILITYARRAY_H_ */
