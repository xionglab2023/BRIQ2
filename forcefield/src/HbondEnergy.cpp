/*
 * HbondEnergy.cpp
 *
 *  Created on: 2023Äê8ÔÂ11ÈÕ
 *      Author: nuc
 */

#include <forcefield/HbondEnergy.h>

namespace NSPforcefield {

HbondEnergy::HbondEnergy() {
	// TODO Auto-generated constructor stub


	this->atLib = new AtomLib();
	ifstream file;
	string fileName;

	fileName = NSPdataio::datapath()+"hbond/hbond.dist";
	file.open(fileName, ios::in);
	if(!file.is_open()) {
		cout << "can't open file: " << fileName << endl;
	}

	int idA, idB;
	int spIndex;

	double dist, logDensityDonor, logDensityAcceptor;
	while(file >> idA >> idB >> dist){
		this->hbDist[idA][idB] = dist;
	}
	file.close();



	fileName = NSPdataio::datapath()+"hbond/relativeHbondLogDensity.txt";
	file.open(fileName, ios::in);
	if(!file.is_open()) {
		cout << "can't open file: " << fileName << endl;
	}

	while(file >> idA >> idB >> spIndex >> logDensityDonor >> logDensityAcceptor){
		this->donorSphereLogDensity[idA][idB][spIndex] = logDensityDonor;
		this->acceptorSphereLogDensity[idA][idB][spIndex] = logDensityAcceptor;
	}
	file.close();

	fileName = NSPdataio::datapath()+"hbond/hbond.ene";
	file.open(fileName, ios::in);
	if(!file.is_open()) {
		cout << "can't open file: " << fileName << endl;
	}

	double ene;
	while(file >> idA >> idB >> ene){
		this->hbEne[idA][idB] = ene;
	}
	file.close();

	fileName = NSPdataio::datapath() + "sphere/sphere5000.index";
	file.open(fileName, ios::in);
	if(!file.is_open()) {
		cout << "can't open file: " << fileName << endl;
	}
	int keyIndex;

	while(file >> keyIndex >> spIndex){
		this->sphereKeyMap[keyIndex] = spIndex;
	}
	file.close();

	fileName = NSPdataio::datapath() + "sphere/sphere5000.index_nearest3Pts";
	file.open(fileName, ios::in);
	if(!file.is_open()) {
		cout << "can't open file: " << fileName << endl;
	}

	int sp1, sp2, sp3;
	double wt1, wt2, wt3;
	int lineId = 0;
	while(file >> keyIndex >> sp1 >> sp2 >> sp3 >> wt1 >> wt2 >> wt3){
		this->keyIdToListID[keyIndex] = lineId;
		this->spIndex1[lineId] = sp1;

		this->spIndex2[lineId] = sp2;
		this->spIndex3[lineId] = sp3;
		this->spWt1[lineId] = wt1;
		this->spWt2[lineId] = wt2;
		this->spWt3[lineId] = wt3;
		lineId ++;
	}
}


double HbondEnergy::getEnergy(int uniqueA, LocalFrame& csA, int uniqueB, LocalFrame& csB){
	int donorIndexA = atLib->uniqueIDToDonorID[uniqueA];
	int acceptorIndexA = atLib->uniqueIDToAcceptorID[uniqueA];

	int donorIndexB = atLib->uniqueIDToDonorID[uniqueB];
	int acceptorIndexB = atLib->uniqueIDToAcceptorID[uniqueB];
	int donorIndex, acceptorIndex;
	LocalFrame csDonor, csAcceptor;

	if(donorIndexA >= 0 && acceptorIndexB >= 0){
		donorIndex = donorIndexA;
		acceptorIndex = acceptorIndexB;
		csDonor = csA;
		csAcceptor = csB;
	}
	else if(donorIndexB >= 0 && acceptorIndexA >= 0){
		donorIndex = donorIndexB;
		acceptorIndex = acceptorIndexA;
		csDonor = csB;
		csAcceptor = csA;
	}
	else {
		return 0.0;
	}

	double e;
	double d0 = hbDist[donorIndex][acceptorIndex];
	double wd = hbEne[donorIndex][acceptorIndex];
	double len = csDonor.origin_.distance(csAcceptor.origin_);
	if(len < d0)
	{
		double u = (len-d0)/d0/0.1;
		e = u*u+wd;
	}
	else
	{
		double u = (len-d0)/d0/0.15;
		if(u>4)
			e = 0.0;
		else
			e = wd*exp(-u*u);
	}


	if(e >= 0.0)
		return e;



	//csDonor: A
	//csAcceptor: B

	double xLen = 1.0/len;

	XYZ tBInA = global2local(csDonor, csAcceptor.origin_)*xLen;
	XYZ tAInB = global2local(csAcceptor, csDonor.origin_)*(xLen);

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

	it = sphereKeyMap.find(idX1*160000+idY1*400+idZ1);
	if(it == sphereKeyMap.end()){
		cout << "len: " << len << endl;
		cout << "B in A: " << tBInA.toString() << endl;
		cout << "invalid sphere index: " << idX1 << " " << idY1 << " " << idZ1 << endl;
		exit(1);
	}
	//int spIndexA = it->second;
	double logDensityA = donorSphereLogDensity[donorIndex][acceptorIndex][it->second];



	it = sphereKeyMap.find(idX2*160000+idY2*400+idZ2);
	if(it == sphereKeyMap.end()){
		cout << "len: " << len << endl;
		cout << "A in B: " << tAInB.toString() << endl;
		cout << "invalid sphere index: " << idX2 << " " << idY2 << " " << idZ2 << endl;
		exit(1);
	}
	//int spIndexB = it->second;
	double logDensityB = acceptorSphereLogDensity[donorIndex][acceptorIndex][it->second];

	//log(1000) = 6.9077
	//need to be verified
    double kOrientation = (logDensityA+logDensityB+log(1000.0))/log(1000.0);

  //  printf("logDenA: %7.3f logDenB: %7.3f kOrientation: %6.4f\n", logDensityA, logDensityB, kOrientation);

    if(kOrientation < 0) return 0.0;
    return e*kOrientation;

}

HbondEnergy::~HbondEnergy() {
	// TODO Auto-generated destructor stub
	delete atLib;
}

} /* namespace NSPforcefield */
