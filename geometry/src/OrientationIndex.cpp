/*
 * OrientationIndex.cpp
 *
 *  Created on: 2023,12,19
 *      Author: pengx
 */

#include <geometry/OrientationIndex.h>

namespace NSPgeometry {

OrientationIndex::OrientationIndex() {
	ifstream file;
	int keyIndex;
	int spIndex;
	double x,y,z;
	string fileName = NSPdataio::datapath() + "sphere/sphere1000.index";

	file.open(fileName, ios::in);
	if(!file.is_open()) {
		cout << "can't open file: " << fileName << endl;
	}

	while(file >> keyIndex >> spIndex){
		this->sphereKeyMap1000[keyIndex] = spIndex;
	}
	file.close();


	fileName = NSPdataio::datapath() + "sphere/sphere2000.index";
	file.open(fileName, ios::in);
	if(!file.is_open()) {
		cout << "can't open file: " << fileName << endl;
	}
	while(file >> keyIndex >> spIndex){
		this->sphereKeyMap2000[keyIndex] = spIndex;
	}
	file.close();

	fileName = NSPdataio::datapath() + "sphere/sphere500.index";
	file.open(fileName, ios::in);
	if(!file.is_open()) {
		cout << "can't open file: " << fileName << endl;
	}
	while(file >> keyIndex >> spIndex){
		this->sphereKeyMap500[keyIndex] = spIndex;
	}
	file.close();

	fileName = NSPdataio::datapath() + "sphere/sphere1000";
	file.open(fileName, ios::in);
	if(!file.is_open()) {
		cout << "can't open file: " << fileName << endl;
	}
	while(file >> x >> y >> z){
		XYZ t(x,y,z);
		this->tList1000.push_back(t);
	}
	file.close();

	fileName = NSPdataio::datapath() + "sphere/sphere2000";
	file.open(fileName, ios::in);
	if(!file.is_open()) {
		cout << "can't open file: " << fileName << endl;
	}
	while(file >> x >> y >> z){
		XYZ t(x,y,z);
		this->tList2000.push_back(t);
	}
	file.close();


	fileName = NSPdataio::datapath() + "sphere/sphere500";
	file.open(fileName, ios::in);
	if(!file.is_open()) {
		cout << "can't open file: " << fileName << endl;
	}
	while(file >> x >> y >> z){
		XYZ t(x,y,z);
		this->tList500.push_back(t);
	}
	file.close();
}


CsMove OrientationIndex::index500ToCsMove(int index){
	int idDist = index/10000000;
	index = index%10000000;
	int idAng = index/250000;
	index = index%250000;
	int spIndexA = index/500;
	int spIndexB = index%500;

	double dist = idDist*0.3;
	double ang = idAng*9;

	LocalFrame csA = getCsA(tList500[spIndexA], ang, dist);
	LocalFrame csB = getCsB(tList500[spIndexB], ang, dist);

	CsMove cm = csB - csA;
	return cm;
}


CsMove OrientationIndex::index1000ToCsMove(int index){
	int idDist = index/40000000;
	index = index%40000000;
	int idAng = index/1000000;
	index = index%1000000;
	int spIndexA = index/1000;
	int spIndexB = index%1000;

	double dist = idDist*0.3;
	double ang = idAng*9;
	LocalFrame csA = getCsA(tList1000[spIndexA], ang, dist);
	LocalFrame csB = getCsB(tList1000[spIndexB], ang, dist);

	CsMove cm = csB - csA;
	return cm;
}

CsMove OrientationIndex::index2000ToCsMove(int indexA, int indexB){
	int idDist = indexA/45;
	int idAng = indexA%45;

	int spIndexA = indexB/2000;
	int spIndexB = indexB%2000;

	double dist = idDist*0.3;
	double ang = idAng*8;

	LocalFrame csA = getCsA(tList2000[spIndexA], ang, dist);
	LocalFrame csB = getCsB(tList2000[spIndexB], ang, dist);

	CsMove cm = csB - csA;
	return cm;
}

int OrientationIndex::moveToIndex500(CsMove& move){
	LocalFrame csA;
	LocalFrame csB = csA + move;

	double len = move.oriMove.length();
	double xLen;
	XYZ tBInA;
	XYZ tAInB;

	if(len == 0) {
		len = 1.0;
		xLen = 1.0;
		tBInA = XYZ(0.0, 0.0, 1.0);
		tAInB = XYZ(0.0, 0.0, -1.0);
	}
	else {
		xLen = 1.0/len;

		tBInA = global2local(csA, csB.origin_)*xLen;
		tAInB = global2local(csB, csA.origin_)*(-xLen);

	}

	int idX1 = (int)((tBInA.x_ + 1.0)*200);
	if(idX1 == 400) idX1 = 399;
	int idY1 = (int)((tBInA.y_ + 1.0)*200);
	if(idY1 == 400) idY1 = 399;
	int idZ1 = (int)((tBInA.z_ + 1.0)*200);
	if(idZ1 == 400) idZ1 = 399;

	int idX2 = (int)((tAInB.x_ + 1.0)*200);
	if(idX2 == 400) idX2 = 399;
	int idY2 = (int)((tAInB.y_ + 1.0)*200);
	if(idY2 == 400) idY2 = 399;
	int idZ2 = (int)((tAInB.z_ + 1.0)*200);
	if(idZ2 == 400) idZ2 = 399;

	map<int, int>::iterator it = sphereKeyMap500.find(idX1*160000+idY1*400+idZ1);
	if(it == sphereKeyMap500.end()){
		cout << "len: " << len << endl;
		cout << "B in A: " << tBInA.toString() << endl;
		cout << "invalid sphere index: " << idX1 << " " << idY1 << " " << idZ1 << endl;
		exit(1);
	}
	int spIndexA = it->second;

	it = sphereKeyMap500.find(idX2*160000+idY2*400+idZ2);
	if(it == sphereKeyMap500.end()){
		cout << "len: " << len << endl;
		cout << "A in B: " << tAInB.toString() << endl;
		cout << "invalid sphere index: " << idX2 << " " << idY2 << " " << idZ2 << endl;
		exit(1);
	}
	int spIndexB = it->second;

	XYZ a = csA.x * (-1.0);
	XYZ b = (csB.origin_ - csA.origin_)*xLen;
	XYZ c = csB.x;
	XYZ n1 = a^b;
	XYZ n2 = b^c;
	XYZ m = n1^b;

	double ang = atan2(m*n2, n1*n2);
	//double ang = 0;
	if(ang < 0)
		ang = -ang * 57.29577952;
	else
		ang = 360 - ang * 57.29577952;

	int indexLen = (int)(len*3.3333333 + 0.5);
	if(indexLen >= 50) indexLen = 49;
	int indexAng = (int)(ang*0.11111111 + 0.5);
	if(indexAng == 40) indexAng = 0;

	return indexLen*10000000+indexAng*250000+spIndexA*500+spIndexB;
}

int OrientationIndex::moveToIndex1000(CsMove& move){
	LocalFrame csA;
	LocalFrame csB = csA + move;

	double len = move.oriMove.length();
	double xLen;
	XYZ tBInA;
	XYZ tAInB;

	if(len == 0) {
		len = 1.0;
		xLen = 1.0;
		tBInA = XYZ(0.0, 0.0, 1.0);
		tAInB = XYZ(0.0, 0.0, -1.0);
	}
	else {
		xLen = 1.0/len;

		tBInA = global2local(csA, csB.origin_)*xLen;
		tAInB = global2local(csB, csA.origin_)*(-xLen);

	}

	int idX1 = (int)((tBInA.x_ + 1.0)*200);
	if(idX1 == 400) idX1 = 399;
	int idY1 = (int)((tBInA.y_ + 1.0)*200);
	if(idY1 == 400) idY1 = 399;
	int idZ1 = (int)((tBInA.z_ + 1.0)*200);
	if(idZ1 == 400) idZ1 = 399;

	int idX2 = (int)((tAInB.x_ + 1.0)*200);
	if(idX2 == 400) idX2 = 399;
	int idY2 = (int)((tAInB.y_ + 1.0)*200);
	if(idY2 == 400) idY2 = 399;
	int idZ2 = (int)((tAInB.z_ + 1.0)*200);
	if(idZ2 == 400) idZ2 = 399;

	map<int, int>::iterator it = sphereKeyMap1000.find(idX1*160000+idY1*400+idZ1);
	if(it == sphereKeyMap1000.end()){
		cout << "len: " << len << endl;
		cout << "B in A: " << tBInA.toString() << endl;
		cout << "invalid sphere index: " << idX1 << " " << idY1 << " " << idZ1 << endl;
		exit(1);
	}
	int spIndexA = it->second;

	it = sphereKeyMap1000.find(idX2*160000+idY2*400+idZ2);
	if(it == sphereKeyMap1000.end()){
		cout << "len: " << len << endl;
		cout << "A in B: " << tAInB.toString() << endl;
		cout << "invalid sphere index: " << idX2 << " " << idY2 << " " << idZ2 << endl;
		exit(1);
	}
	int spIndexB = it->second;

	XYZ a = csA.x * (-1.0);
	XYZ b = (csB.origin_ - csA.origin_)*xLen;
	XYZ c = csB.x;
	XYZ n1 = a^b;
	XYZ n2 = b^c;
	XYZ m = n1^b;

	double ang = atan2(m*n2, n1*n2);
	//double ang = 0;
	if(ang < 0)
		ang = -ang * 57.29577952;
	else
		ang = 360 - ang * 57.29577952;

	int indexLen = (int)(len*3.3333333 + 0.5);
	if(indexLen >= 50) indexLen = 49;
	int indexAng = (int)(ang*0.11111111 + 0.5);
	if(indexAng == 40) indexAng = 0;

	return indexLen*40000000+indexAng*1000000+spIndexA*1000+spIndexB;
}

pair<int,int> OrientationIndex::moveToIndex2000(CsMove& move){
	LocalFrame csA;
	LocalFrame csB = csA + move;

	double len = move.oriMove.length();
	double xLen;
	XYZ tBInA;
	XYZ tAInB;

	if(len == 0) {
		len = 1.0;
		xLen = 1.0;
		tBInA = XYZ(0.0, 0.0, 1.0);
		tAInB = XYZ(0.0, 0.0, -1.0);
	}
	else {
		xLen = 1.0/len;

		tBInA = global2local(csA, csB.origin_)*xLen;
		tAInB = global2local(csB, csA.origin_)*(-xLen);

	}

	int idX1 = (int)((tBInA.x_ + 1.0)*200);
	if(idX1 == 400) idX1 = 399;
	int idY1 = (int)((tBInA.y_ + 1.0)*200);
	if(idY1 == 400) idY1 = 399;
	int idZ1 = (int)((tBInA.z_ + 1.0)*200);
	if(idZ1 == 400) idZ1 = 399;

	int idX2 = (int)((tAInB.x_ + 1.0)*200);
	if(idX2 == 400) idX2 = 399;
	int idY2 = (int)((tAInB.y_ + 1.0)*200);
	if(idY2 == 400) idY2 = 399;
	int idZ2 = (int)((tAInB.z_ + 1.0)*200);
	if(idZ2 == 400) idZ2 = 399;

	map<int, int>::iterator it = sphereKeyMap2000.find(idX1*160000+idY1*400+idZ1);
	if(it == sphereKeyMap2000.end()){
		cout << "len: " << len << endl;
		cout << "B in A: " << tBInA.toString() << endl;
		cout << "invalid sphere index: " << idX1 << " " << idY1 << " " << idZ1 << endl;
		exit(1);
	}
	int spIndexA = it->second;

	it = sphereKeyMap2000.find(idX2*160000+idY2*400+idZ2);
	if(it == sphereKeyMap2000.end()){
		cout << "len: " << len << endl;
		cout << "A in B: " << tAInB.toString() << endl;
		cout << "invalid sphere index: " << idX2 << " " << idY2 << " " << idZ2 << endl;
		exit(1);
	}
	int spIndexB = it->second;

	XYZ a = csA.x * (-1.0);
	XYZ b = (csB.origin_ - csA.origin_)*xLen;
	XYZ c = csB.x;
	XYZ n1 = a^b;
	XYZ n2 = b^c;
	XYZ m = n1^b;

	double ang = atan2(m*n2, n1*n2);
	//double ang = 0;
	if(ang < 0)
		ang = -ang * 57.29577952;
	else
		ang = 360 - ang * 57.29577952;

	int indexLen = (int)(len*3.3333333 + 0.5);
	if(indexLen >= 50) indexLen = 49;
	int indexAng = (int)(ang*0.125 + 0.5);
	if(indexAng == 45) indexAng = 0;

	int indexA = indexLen*45+indexAng;
	int indexB = spIndexA*2000+spIndexB;
	pair<int, int> p(indexA, indexB);
	return p;
}

OrientationIndex::~OrientationIndex() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPpredNA */
