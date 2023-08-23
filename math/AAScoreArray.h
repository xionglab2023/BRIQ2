/*
 * AAScoreArray.h
 *
 *  Created on: 2021Äê12ÔÂ10ÈÕ
 *      Author: pengx
 */

#ifndef MATH_AASCOREARRAY_H_
#define MATH_AASCOREARRAY_H_

#include <stdlib.h>
#include <string>
#include <iostream>
#include <math.h>
#include "tools/StringTool.h"
#include "math/AAProbabilityArray.h"

namespace NSPmath {

using namespace std;

class AAScoreArray {

public:
	double sa[20];
	double saBg[20];

	AAScoreArray();
	AAScoreArray(const string& line);
	AAScoreArray(AAProbabilityArray* pa);

	AAScoreArray& operator=(const AAScoreArray& other){
		for(int i=0;i<20;i++){
			this->sa[i] = other.sa[i];
		}
		return *this;
	}

	void multipy(float wt);

	void print(ofstream& out){
		for(int i=0;i<19;i++){
			out << sa[i] << " ";
		}
		out << sa[19] << endl;
	}

	virtual ~AAScoreArray();
};

} /* namespace NSPmath */

#endif /* MATH_AASCOREARRAY_H_ */
