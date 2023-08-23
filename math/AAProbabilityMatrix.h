/*
 * AAProbabilityMatrix.h
 *
 */

#ifndef MATH_AAPROBABILITYMATRIX_H_
#define MATH_AAPROBABILITYMATRIX_H_

#include <stdio.h>
#include <iostream>
namespace NSPmath {

using namespace std;

class AAProbabilityMatrix {
public:
	double pm[20][20];
	double sampleNum;

	AAProbabilityMatrix();
	AAProbabilityMatrix(double *p, double N);
	AAProbabilityMatrix(int *c);
	void normalize();
	void addPseudoCount(AAProbabilityMatrix* bg, int pseudoNum);
	void print();

	virtual ~AAProbabilityMatrix();
};

} /* namespace NSPmath */

#endif /* MATH_AAPROBABILITYMATRIX_H_ */
