/*
 * AAProbabilityArray.cpp
 */

#include "math/AAProbabilityArray.h"

namespace NSPmath {

AAProbabilityArray::AAProbabilityArray() {
	// TODO Auto-generated constructor stub
	for(int i=0;i<20;i++){
		this->pa[i] = 0.05;
	}
}

AAProbabilityArray::AAProbabilityArray(double* p, double N){
	this->sampelNum = N;
	for(int i=0;i<20;i++){
		this->pa[i] = p[i];
	}
	normalize();
}

AAProbabilityArray::AAProbabilityArray(int* count){
	this->sampelNum = 0;
	for(int i=0;i<20;i++){
		this->pa[i] = count[i];
		this->sampelNum += count[i];
	}
	normalize();
}

void AAProbabilityArray::normalize(){
	double sum = 0;
	for(int i=0;i<20;i++){
		sum += pa[i];
	}
	if(sum != 0) {
		for(int i=0;i<20;i++){
			pa[i] = pa[i]/sum;
		}
	}
}

void AAProbabilityArray::addPseudoCount(SubstitutionMatrix* subM, int pseudoNum){
	double pb[20];
	double sum = 0;
	for(int i=0;i<20;i++){
		pb[i] = 0.0;
	}
	for(int i=0;i<20;i++){
		for(int j=0;j<20;j++){
			pb[i]+=pa[j]*subM->subM[i][j];
		}
		sum += pb[i];
	}

	double n = this->sampelNum - 1;
	if(n < 0)
		n = 20;

	for(int i=0;i<20;i++){
		this->pa[i] = (n*pa[i] + pseudoNum*pb[i])/(n+pseudoNum);
	}

	normalize();
}

double AAProbabilityArray::entropy(){
	double ent = 0;
	double log2 = log(2.0);
	for(int i=0;i<20;i++){
		if(pa[i] > 0){
			ent -= pa[i]*log(pa[i])/log2;
		}
	}
	return ent;
}

double AAProbabilityArray::distance(AAProbabilityArray* other){
	double D = 0;
	double log2 = log(2.0);
	for(int i=0;i<20;i++)
	{
		double p=pa[i];
		double q=other->pa[i];
		if(p!=0)
		{
			D+=p*log(2*p/(p+q))/log2;
		}

		if(q!=0)
		{
			D+=q*log(2*q/(p+q))/log2;
		}
	}
	return D;
}

void AAProbabilityArray::printAAName(){
	NSPmodel::ResName rn;
	for(int i=0;i<20;i++)
	{
		printf("  %s  ", rn.intToTri(i).c_str());
	}
	printf("\n");
}

void AAProbabilityArray::print(){
	for(int i=0;i<20;i++){
		printf("%6.4f ", pa[i]);
	}
	printf("\n");
}

AAProbabilityArray::~AAProbabilityArray() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPmath */
