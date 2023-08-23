/*
 * TestTime.cpp
 *
 */

#include <time.h>
#include <stdlib.h>
#include <iostream>
#include "geometry/xyz.h"
#include "geometry/CsMove.h"
#include "geometry/localframe.h"

using namespace std;
using namespace NSPgeometry;

int main(){
	clock_t start = clock();


	cout << "start" << endl;
	int x[30000][100];
	for(int i=0;i<30000;i++) {
		for(int j=0;j<100;j++) {
			x[i][j] = 0;
		}
	}
	cout << "finish init" << endl;


	XYZ a(0, 0, 0);
	XYZ b(2, 5, 9);
	XYZ c(6, 0, 8);

	XYZ d(2, 1, 3);
	LocalFrame csA(a, b, c);
	LocalFrame csB(b,c,d);
	CsMove mv = csB - csA;
	//CsMove mv2;
	for(int i=0;i<100000000;++i) {
		//XYZ e = local2global(csA, d);
		//LocalFrame lf = csA + mv;
		LocalFrame lf = csA.add(mv);
		//double x = a-b;
		//double y = c-d;
		//double z = e-f;
		//double dd = x*x + y*y + z*z;
	}

	clock_t end = clock();
	cout << "time: " << (float)(end-start)/CLOCKS_PER_SEC << "s" << endl;
}

