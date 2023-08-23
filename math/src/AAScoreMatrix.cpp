/*
 * AAScoreMatrix.cpp
 *
 */

#include "math/AAScoreMatrix.h"

namespace NSPmath {

AAScoreMatrix::AAScoreMatrix() {
	for(int i=0;i<20;i++){
		for(int j=0;j<20;j++){
			sm[i][j] = 0.0;
		}
	}
	this->sampleNum = -1;
}

AAScoreMatrix::AAScoreMatrix(double* sm, double N){
	for(int i=0;i<20;i++){
		for(int j=0;j<20;j++){
			this->sm[i][j] = sm[i*20+j];
		}
	}
	this->sampleNum = N;
}

AAScoreMatrix::AAScoreMatrix(AAProbabilityMatrix* pm, AAProbabilityMatrix* pmBg){
	for(int i=0;i<20;i++){
		for(int j=0;j<20;j++){
			double p = pm->pm[i][j];
			double bg = pmBg->pm[i][j];
			if(p == 0 && bg == 0)
				sm[i][j] = 0.0;
			else if(p==0 && bg >0)
				sm[i][j] = 20.0;
			else if(p>0 && bg == 0)
				sm[i][j] = -20.0;
			else {
				double s = log(bg/p);
				if(s > 20)
					s = 20;
				if(s < -20)
					s = -20;
				sm[i][j] = s;
			}
		}
	}
	this->sampleNum = pm->sampleNum;
}

void AAScoreMatrix::multipy(double wt){
	for(int i=0;i<20;i++){
		for(int j=0;j<20;j++){
			this->sm[i][j] = this->sm[i][j]*wt;
		}
	}
}

void AAScoreMatrix::print(ofstream& out){
	char xx[200];
	sprintf(xx, "%-7.1f", this->sampleNum);
	out<< string(xx) << endl;
	for(int i=0;i<20;i++){
		for(int j=0;j<20;j++){
			sprintf(xx, "%-7.4f ", sm[i][j]);
			out<< string(xx);
		}
		out<< endl;
	}
}

AAScoreMatrix::~AAScoreMatrix() {
}

} /* namespace NSPmath */
