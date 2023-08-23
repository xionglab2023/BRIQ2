/*
 * TestKey.cpp
 *
 */




#include "model/BaseRotamer.h"
#include "model/BaseRotamerLib.h"
#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>

#include "model/StructureModel.h"
#include "predNA/MCRun.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredna;
using namespace std;

LocalFrame getCsA(XYZ t, double dihed, double dist){
	double x = t.x_;
	double y = t.y_;
	double z;
	z = sqrt(1-x*x-y*y);
	if(t.z_ < 0) z = -z;
	double sx = sqrt(1-x*x);
	double sy = sqrt(1-y*y);
	double sz = sqrt(1-z*z);

	double sinD1, cosD1, sinD2, cosD2;

	sinD1 = z/sx/sy;
	cosD1 = -x*y/sx/sy;


	double ang1 = atan2(sinD1, cosD1);

	double ang0 = dihed*0.008726646;
	double tm[3][3];
	tm[0][0] = x;
	tm[0][1] = y;
	tm[0][2] = z;
	tm[1][0] = sx*sin(ang0);
	tm[1][1] = sy*sin(ang0-ang1);
	tm[2][0] = sx*cos(ang0);
	tm[2][1] = sy*cos(ang0-ang1);

	XYZ t1(x, tm[1][0], tm[2][0]);
	XYZ t2(y, tm[1][1], tm[2][1]);
	XYZ t3 = t1^t2;
	tm[1][2] = t3.y_;
	tm[2][2] = t3.z_;

	XYZ ori(-dist*0.5, 0, 0);
	LocalFrame cs(ori, tm);
	return cs;
}

LocalFrame getCsB(XYZ t, double dihed, double dist){
	double x = t.x_;
	double y = t.y_;
	double z;
	z = sqrt(1-x*x-y*y);
	if(t.z_ < 0) z = -z;
	double sx = sqrt(1-x*x);
	double sy = sqrt(1-y*y);
	double sz = sqrt(1-z*z);

	double sinD1, cosD1, sinD2, cosD2;

	sinD1 = z/sx/sy;
	cosD1 = -x*y/sx/sy;


	double ang1 = atan2(sinD1, cosD1);

	double ang0 = -dihed*0.008726646;
	double tm[3][3];
	tm[0][0] = x;
	tm[0][1] = y;
	tm[0][2] = z;
	tm[1][0] = sx*sin(ang0);
	tm[1][1] = sy*sin(ang0-ang1);
	tm[2][0] = sx*cos(ang0);
	tm[2][1] = sy*cos(ang0-ang1);

	XYZ t1(x, tm[1][0], tm[2][0]);
	XYZ t2(y, tm[1][1], tm[2][1]);
	XYZ t3 = t1^t2;
	tm[1][2] = t3.y_;
	tm[2][2] = t3.z_;

	XYZ ori(dist*0.5, 0, 0);
	LocalFrame cs(ori, tm);
	return cs;
}


int main(int argc, char** argv){




	ifstream file;

	file.open("/export/home/s2982206/aspen/cppWorkplace/RNAModeling/data/sphere800", ios::in);

	vector<XYZ> tList;
	double x,y,z;
	while(file >> x >> y >> z){
		XYZ t(x,y,z);
		tList.push_back(t);
	}
	file.close();

	char ss[7];


	int i=100;
	int j = 200;
	int k = 20;
	int l = 30;

	ss[0] = i/30 + '!';
	ss[1] = i%30 + '!';
	ss[2] = j/30 + '!';
	ss[3] = j%30 + '!';
	ss[4] = k + '!';
	ss[5] = l + '!';

	XYZ t1 = tList[i];
	XYZ t2 = tList[j];
	double d = k*0.25;
	double ang = l*7.2;

	LocalFrame csA = getCsA(t1, ang, d);
	LocalFrame csB = getCsB(t2, ang, d);
	string key = string(ss);
	CsMove mv = csB - csA;
	printf("%5.3f %8.3f\n", d, ang);
	cout << key << " " << mv.toString() << endl;

}

