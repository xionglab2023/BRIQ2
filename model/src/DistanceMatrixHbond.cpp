/*
 * DistanceMatrixHbond.cpp
 *
 */

#include <model/DistanceMatrixHbond.h>

namespace NSPmodel {

DistanceMatrixHbond::DistanceMatrixHbond() {
	// TODO Auto-generated constructor stub
	for(int i=0;i<16;i++){
		this->dm[i] = 0;
	}
}

DistanceMatrixHbond::DistanceMatrixHbond(const string& line){
	vector<string> list;
	NSPtools::splitString(line, " ", &list);
	double d;
	for(int i=0;i<16;i++){
		d = atof(list[i+2].c_str());
		this->dm[i] = d;
	}
}

DistanceMatrixHbond::DistanceMatrixHbond(LocalFrame& csA, LocalFrame& csB){
	XYZ a(3.464 , 0 ,  -1.224);
	XYZ b(-1.732 ,  3.0 , -1.224);
	XYZ c(-1.732 ,  -3.0 , -1.224);
	XYZ d(0,  0,  3.674);

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

double DistanceMatrixHbond::distanceTo(DistanceMatrixHbond& other){
	double dd = 0;
	for(int i=0;i<16;i++){
		dd += (this->dm[i] - other.dm[i])* (this->dm[i] - other.dm[i]);
	}
	return sqrt(dd*0.0625);
}

double DistanceMatrixHbond::squareDistance(DistanceMatrixHbond& other){
	double dd = 0;
	for(int i=0;i<16;i++){
		dd += (this->dm[i] - other.dm[i])* (this->dm[i] - other.dm[i]);
	}
	return dd*0.0625;
}


DistanceMatrixHbond::~DistanceMatrixHbond() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPmodel */
