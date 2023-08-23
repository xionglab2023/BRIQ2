/*
 * TestDM.cpp
 *
 *  Created on: Apr 16, 2019
 *      Author: s2982206
 */


#include "model/BaseDistanceMatrix.h"

using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;

int main(){

	XYZ a(0, 0, 0);
	XYZ b(2, 5, 9);
	XYZ c(6, 0, 8);

	XYZ d(-5, -2, -2);
	LocalFrame csA(a, b, c);
	LocalFrame csB(d, c, b);
	BaseDistanceMatrix dm(csA, csB);
	cout << dm.hashKey << endl;

}


