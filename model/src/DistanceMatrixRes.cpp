/*
 * DistanceMatrixRes.cpp
 *
 */

#include "model/DistanceMatrixRes.h"

namespace NSPmodel {

DistanceMatrixRes::DistanceMatrixRes() {
	// TODO Auto-generated constructor stub
	for(int i=0;i<16;i++){
		this->dm[i] = 0;
	}
}

DistanceMatrixRes::DistanceMatrixRes(const string& line){
	vector<string> list;
	NSPtools::splitString(line, " ", &list);
	double d;
	for(int i=0;i<16;i++){
		d = atof(list[i+2].c_str());
		this->dm[i] = d;
	}
}

DistanceMatrixRes::DistanceMatrixRes(LocalFrame& csA, LocalFrame& csB){
	XYZ a( 2.397, -0.218, -0.449);
	XYZ b(-2.799,  2.782, -0.449);
	XYZ c(-2.799, -3.218, -0.449);
	XYZ d(-1.067, -0.218,  4.449);

	XYZ a1 = local2global(csA, a);
	XYZ b1 = local2global(csA, b);
	XYZ c1 = local2global(csA, c);
	XYZ d1 = local2global(csA, d);

	XYZ a2 = local2global(csB, a);
	XYZ b2 = local2global(csB, b);
	XYZ c2 = local2global(csB, c);
	XYZ d2 = local2global(csB, d);

	dm[0] = a1.distance(a2);
	dm[1] = a1.distance(b2);
	dm[2] = a1.distance(c2);
	dm[3] = a1.distance(d2);
	dm[4] = b1.distance(a2);
	dm[5] = b1.distance(b2);
	dm[6] = b1.distance(c2);
	dm[7] = b1.distance(d2);
	dm[8] = c1.distance(a2);
	dm[9] = c1.distance(b2);
	dm[10] = c1.distance(c2);
	dm[11] = c1.distance(d2);
	dm[12] = d1.distance(a2);
	dm[13] = d1.distance(b2);
	dm[14] = d1.distance(c2);
	dm[15] = d1.distance(d2);
}


DistanceMatrixRes::DistanceMatrixRes(CsMove& cm){
	LocalFrame csA;
	LocalFrame csB = csA.add(cm);
	XYZ a( 2.397, -0.218, -0.449);
	XYZ b(-2.799,  2.782, -0.449);
	XYZ c(-2.799, -3.218, -0.449);
	XYZ d(-1.067, -0.218,  4.449);

	XYZ a2 = local2global(csB, a);
	XYZ b2 = local2global(csB, b);
	XYZ c2 = local2global(csB, c);
	XYZ d2 = local2global(csB, d);

	dm[0] = a.distance(a2);
	dm[1] = a.distance(b2);
	dm[2] = a.distance(c2);
	dm[3] = a.distance(d2);
	dm[4] = b.distance(a2);
	dm[5] = b.distance(b2);
	dm[6] = b.distance(c2);
	dm[7] = b.distance(d2);
	dm[8] = c.distance(a2);
	dm[9] = c.distance(b2);
	dm[10] = c.distance(c2);
	dm[11] = c.distance(d2);
	dm[12] = d.distance(a2);
	dm[13] = d.distance(b2);
	dm[14] = d.distance(c2);
	dm[15] = d.distance(d2);
}

double DistanceMatrixRes::distanceTo(DistanceMatrixRes& other){
	double dd = 0;
	for(int i=0;i<16;i++){
		dd += (this->dm[i] - other.dm[i])* (this->dm[i] - other.dm[i]);
	}
	return sqrt(dd*0.0625);
}

double DistanceMatrixRes::squareDistance(DistanceMatrixRes& other){
	double dd = 0;
	for(int i=0;i<16;i++){
		dd += (this->dm[i] - other.dm[i])* (this->dm[i] - other.dm[i]);
	}
	return dd*0.0625;
}

DistanceMatrixRes::~DistanceMatrixRes() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPmodel */
