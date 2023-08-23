/*
 * AAScoreMatrix.h
 *
 */

#ifndef MATH_AASCOREMATRIX_H_
#define MATH_AASCOREMATRIX_H_

#include <string>
#include <math.h>
#include "math/AAProbabilityMatrix.h"
#include <iostream>
#include <fstream>

namespace NSPmath {

using namespace std;

class AAScoreMatrix {
public:
	double sm[20][20];
	double sampleNum;

	AAScoreMatrix();
	AAScoreMatrix(double *sm, double N);
	AAScoreMatrix(AAProbabilityMatrix* pm, AAProbabilityMatrix* pmBg);
	AAScoreMatrix& operator=(const AAScoreMatrix& other){
		int i,j;
		for(i=0;i<20;i++){
			for(j=0;j<20;j++){
				this->sm[i][j] = other.sm[i][j];
			}
		}
		this->sampleNum = other.sampleNum;
		return *this;
	}

	void multipy(double wt);
	void print(ofstream& out);

	virtual ~AAScoreMatrix();
};

} /* namespace NSPmath */

#endif /* MATH_AASCOREMATRIX_H_ */
