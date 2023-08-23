/*
 * Stat.cpp
 *
 */

#include "math/Stat.h"

namespace NSPmath {

double sd(vector<double>& dList){
	int n = dList.size();
	if(n == 0) return 0.0;

	double sum = 0;
	double squareSum = 0;
	for(int i=0;i<n;i++){
		sum += dList[i];
		squareSum += dList[i]*dList[i];
	}
	double mean = sum/n;
	return sqrt(squareSum/n - mean*mean);
}

double mid(vector<double>& dList){
	int n = dList.size();
	if(n == 0) return 0.0;
	double array[dList.size()];
	for(int i=0;i<n;i++){
		array[i] = dList[i];
	}
	sort(array, array+n);
	return array[n/2];
}

double mean(vector<double>& dList){
	int n = dList.size();
	if(n == 0) return 0.0;

	double sum = 0;
	for(int i=0;i<n;i++){
		sum += dList[i];
	}
	return sum/n;
}

double sum(vector<double>& dList){
	int n = dList.size();
	if(n == 0) return 0.0;
	double sum = 0;
	for(int i=0;i<n;i++){
		sum += dList[i];
	}
	return sum;
}

double top(vector<double>& dList, double p){
	int n = dList.size();
	if(n == 0) return 0.0;

	double array[dList.size()];
	for(int i=0;i<n;i++){
		array[i] = dList[i];
	}
	sort(array, array+n);

	int index = (int)n*p;
	if(index< 0) index = 0;
	if(index >= n) index = n-1;
	return array[index];
}

float top(vector<float>& dList, double p){
	int n = dList.size();
	if(n == 0) return 0.0;

	float array[dList.size()];
	for(int i=0;i<n;i++){
		array[i] = dList[i];
	}
	sort(array, array+n);

	int index = (int)n*p;
	if(index< 0) index = 0;
	if(index >= n) index = n-1;
	return array[index];
}


} /* namespace NSPmath */
