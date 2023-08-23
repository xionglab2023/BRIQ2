/*
 * TestCoordSystem.cpp
 *
 */

#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include "geometry/TransMatrix.h"
#include "geometry/CsMove.h"
#include "geometry/RMSD.h"
#include "geometry/Angles.h"
#include <iostream>

using namespace std;
using namespace NSPgeometry;

LocalFrame getCsA(double x, double y, double z, double dihed, double dist){
	double sx = sqrt(1-x*x);
	double sy = sqrt(1-y*y);
	double sz = sqrt(1-z*z);
	double ang1 = asin(-x*y/sx/sy);
	double ang2 = asin(-x*z/sx/sz);
	double ang0 = dihed*0.008726646;
	double tm[3][3];
	tm[0][0] = x;
	tm[0][1] = y;
	tm[0][2] = z;
	tm[1][0] = sx*sin(ang0);
	tm[1][1] = sy*cos(ang1-ang0);
	tm[1][2] = sz*cos(ang2-ang0);
	tm[2][0] = sx*cos(ang0);
	tm[2][1] = sy*sin(ang1-ang0);
	tm[2][2] = sz*sin(ang2-ang0);

	XYZ ori(-dist*0.5, 0, 0);
	LocalFrame cs(ori, tm);
	return cs;
}

LocalFrame getCsB(double x, double y, double z, double dihed, double dist){
	double sx = sqrt(1-x*x);
	double sy = sqrt(1-y*y);
	double sz = sqrt(1-z*z);
	double ang1 = asin(-x*y/sx/sy);
	double ang2 = asin(-x*z/sx/sz);
	double ang0 = -dihed*0.008726646;
	double tm[3][3];
	tm[0][0] = -x;
	tm[0][1] = -y;
	tm[0][2] = -z;
	tm[1][0] = sx*sin(ang0);
	tm[1][1] = sy*cos(ang1-ang0);
	tm[1][2] = sz*cos(ang2-ang0);
	tm[2][0] = sx*cos(ang0);
	tm[2][1] = sy*sin(ang1-ang0);
	tm[2][2] = sz*sin(ang2-ang0);

	XYZ ori(dist*0.5, 0, 0);
	LocalFrame cs(ori, tm);
	return cs;
}


int main(){
	{

		double x1 = 0.79319;
		double y1 = 0.36201;
		double z1 = -0.48969;

		double x2 = 0.94328;
		double y2 = -0.255;
		double z2 = -0.2126;


		LocalFrame csA = getCsA(x1,y1,z1,-45,6.0);
		LocalFrame csB = getCsB(x2,y2,z2,-45,6.0);

		XYZ oriA = csA.origin_;
		XYZ oriB = csB.origin_;

		XYZ oriAInB = csB.global2localcrd(oriA);
		XYZ oriBInA = csA.global2localcrd(oriB);

		XYZ la = oriBInA/6.0;
		XYZ lb = oriAInB/6.0;

		XYZ t(1.0, 0, 0);
		XYZ t1 = csA.local2globalcrd(t);
		XYZ t2 = csB.local2globalcrd(t);

		double dihed = dihedral(t1, csA.origin_,  csB.origin_, t2);
		cout << la.toString() << endl;
		cout << lb.toString() << endl;
		cout << dihed << endl;

		XYZ A(0,0,0);
		XYZ B(1,1,0);
		XYZ C(0, 2, 0);
		XYZ D(1.4, 2.1, 3.0);
		XYZ E(3.12, 4.11, 9.35);

		LocalFrame cs1(A, B, C);
		LocalFrame cs2(B, C, D);
		LocalFrame cs3(C, D, E);
		CsMove mv1 = cs1 - cs2;
		CsMove mv2 = cs3 - cs2;

		CsMove mv3 = mv1 - mv2;
		CsMove mv4 = mv3 + mv2;
		cout << "move1" << endl;
		mv1.print();
		cout << "move3" << endl;
		mv3.print();
		cout << "move4" << endl;
		mv4.print();

	}


}


