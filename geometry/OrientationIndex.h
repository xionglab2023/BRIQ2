/*
 * OrientationIndex.h
 *
 *  Created on: 2023,12,19
 *      Author: pengx
 */

#ifndef PREDNA_ORIENTATIONINDEX_H_
#define PREDNA_ORIENTATIONINDEX_H_

#include <vector>
#include "geometry/localframe.h"
#include "geometry/xyz.h"
#include "geometry/Angles.h"
#include <fstream>
#include <map>
#include "geometry/CsMove.h"
#include "dataio/datapaths.h"
#include "tools/StringTool.h"


namespace NSPgeometry {

using namespace std;


class OrientationIndex {

private:

	/*
	 * length: 0~15
	 * dihedral angles: -180~180
	 * sphereA:
	 * sphereB:
	 */
	map<int, int> sphereKeyMap500; //50*40*500*500
	map<int, int> sphereKeyMap1000; //50*40*1000*1000
	map<int, int> sphereKeyMap2000; //50*45*2000*2000

	vector<XYZ> tList500;
	vector<XYZ> tList1000;
	vector<XYZ> tList2000;


public:
	OrientationIndex(bool withBinary=false);
	int dump(ostream& outs);
	int load(istream& ins);
	const map<int, int>& getSKM1000() const {
		return sphereKeyMap1000;
	}
	CsMove index500ToCsMove(int index);
	CsMove index1000ToCsMove(int index);
	CsMove index1000ToCsMoveWithRandPerturbation(int index);
	CsMove fixIndex1000WithRandPerturbation(CsMove& cm);
	CsMove index2000ToCsMove(int indexA, int indexB);


	int moveToIndex500(CsMove& move);
	int moveToIndex1000(CsMove& move);
	pair<int, int> moveToIndex2000(CsMove& move);
	virtual ~OrientationIndex();
};


inline LocalFrame getCsA(XYZ t, double dihed, double dist){
	double x = t.x_;
	double y = t.y_;
	double z;
	z = sqrt(1-x*x-y*y);
	if(1-x*x-y*y <= 0)
		z = 0.0000001;
	if(t.z_ < 0) z = -z;
	double sx = sqrt(1-x*x);
	double sy = sqrt(1-y*y);
	double sz = sqrt(1-z*z);

	if(1-x*x <= 0)
		sx = 0.0000001;
	if(1-y*y <= 0)
		sy = 0.0000001;
	if(1-z*z <= 0)
		sz = 0.0000001;

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

inline LocalFrame getCsB(XYZ t, double dihed, double dist){
	double x = t.x_;
	double y = t.y_;
	double z;
	z = sqrt(1-x*x-y*y);
	
	if(1-x*x-y*y <= 0)
		z = 0.0000001;
	
	if(t.z_ < 0) z = -z;
	double sx = sqrt(1-x*x);
	double sy = sqrt(1-y*y);
	double sz = sqrt(1-z*z);

	if(1-x*x <= 0)
		sx = 0.0000001;
	if(1-y*y <= 0)
		sy = 0.0000001;
	if(1-z*z <= 0)
		sz = 0.0000001;

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
} /* namespace NSPpredNA */

#endif /* PREDNA_ORIENTATIONINDEX_H_ */
