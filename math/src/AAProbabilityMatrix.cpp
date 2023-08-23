/*
 * AAProbabilityMatrix.cpp
 *
 */

#include "math/AAProbabilityMatrix.h"

namespace NSPmath {

AAProbabilityMatrix::AAProbabilityMatrix() {
	for(int i=0;i<20;i++){
		for(int j=0;j<20;j++){
			this->pm[i][j] = 0.0025;
		}
	}
	this->sampleNum = -1;
}

AAProbabilityMatrix::AAProbabilityMatrix(double *p, double N){
	for(int i=0;i<20;i++){
		for(int j=0;j<20;j++){
			this->pm[i][j] = p[i*20+j];
		}
	}
	this->sampleNum = N;
	normalize();
}

AAProbabilityMatrix::AAProbabilityMatrix(int *p){
	this->sampleNum = 0;
	for(int i=0;i<20;i++){
		for(int j=0;j<20;j++){
			this->pm[i][j] = p[i*20+j];
			this->sampleNum += 1.0;
		}
	}
	normalize();
}

void AAProbabilityMatrix::normalize(){
	double sum = 0.0;
	for(int i=0;i<20;i++){
		for(int j=0;j<20;j++){
			sum += pm[i][j];
		}
	}
	if(sum == 0) return;
	for(int i=0;i<20;i++){
		for(int j=0;j<20;j++){
			pm[i][j] = pm[i][j]/sum;
		}
	}
}

void AAProbabilityMatrix::addPseudoCount(AAProbabilityMatrix* bg, int pseudoNum){
	if(this->sampleNum < 20){
		for(int i=0;i<20;i++){
			for(int j=0;j<20;j++){
				this->pm[i][j] = bg->pm[i][j];
			}
		}
	}
	else {
		for(int i=0;i<20;i++){
			for(int j=0;j<20;j++){
				this->pm[i][j] = pm[i][j]*(sampleNum-20) + bg->pm[i][j]*pseudoNum;
			}
		}
	}
	normalize();
}


void AAProbabilityMatrix::print(){
	for(int i=0;i<20;i++)
	{
		for(int j=0;j<20;j++)
		{
			printf("%6.4f ", pm[i][j]);
		}
		printf("\n");
	}
}

AAProbabilityMatrix::~AAProbabilityMatrix() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPmath */
