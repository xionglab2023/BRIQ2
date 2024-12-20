
#include <forcefield/Xtb6dEnergy.h>

namespace NSPforcefield {

CsMoveTo6DKeySP500::CsMoveTo6DKeySP500(){
	this->spherePointNum = 500;
	ifstream file;
	char xx[20];
	sprintf(xx, "%d", spherePointNum);
	string fileName = NSPdataio::datapath() + "sphere/sphere" + string(xx) + ".index";
	file.open(fileName, ios::in);
	if(!file.is_open()) {
		cout << "can't open file: " << fileName << endl;
	}
	int keyIndex;
	int spIndex;
	while(file >> keyIndex >> spIndex){
		this->sphereKeyMap[keyIndex] = spIndex;
	}
	file.close();

	fileName = NSPdataio::datapath() + "sphere/sphere" + string(xx) + ".index_nearest3Pts";
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

	int indexA, indexB, indexC, indexD;
	double dA, dB, dC, dD;
	double wA, wB, wC, wD;
	for(int i=0;i<500;i++){
		double d = i*0.03+0.015;
		if(i >= 490){
			indexA = 49;
			indexB = 0;
			dA = (d-indexA*0.3)*3.3333333333;
			dB = 1-dA;
		}
		else {
			indexA = i/10;
			indexB = indexA+1;
			dA = (d-indexA*0.3)*3.3333333333;
			dB = 1-dA;
		}
		for(int j=0;j<400;j++){
			double ang = j*0.9+0.045;
			if(j >= 390){
				indexC = 39;
				indexD = 0;
				dC = (ang - indexC*9)*0.111111;
				dD = 1 - dC;
			}
			else {
				indexC = j/10;
				indexD = indexC+1;
				dC = (ang - indexC*9)*0.111111;
				dD = 1 - dC;
			}

			//indexA, indexB: index for distance
			//indexC, indexD: index for angle

			this->biIndex[i*1600 + j*4] = indexA*40+indexC;
			this->biIndex[i*1600 + j*4 + 1] = indexA*40+indexD;
			this->biIndex[i*1600 + j*4 + 2] = indexB*40+indexC;
			this->biIndex[i*1600 + j*4 + 3] = indexB*40+indexD;

			this->biWt[i*1600 + j*4] = dB * dD;
			this->biWt[i*1600 + j*4 + 1] = dB * dC;
			this->biWt[i*1600 + j*4 + 2] = dA * dD;
			this->biWt[i*1600 + j*4 + 3] = dA * dC;
		}
	}
}

void CsMoveTo6DKeySP500::getIndexAndWeight(const LocalFrame& csA, const LocalFrame& csB, double len, int outIndex[5], double outWeights[4]){

	if(len == 0) {
		int index = 0;

		outIndex[0] = this->biIndex[index];
		outIndex[1] = this->biIndex[index+1];
		outIndex[2] = this->biIndex[index+2];
		outIndex[3] = this->biIndex[index+3];

		outIndex[4] = 0;
		outWeights[0] = this->biWt[index];
		outWeights[1] = this->biWt[index+1];
		outWeights[2] = this->biWt[index+2];
		outWeights[3] = this->biWt[index+3];
		return;
	}

	double xLen = 1.0/len;

	XYZ tBInA = global2local(csA, csB.origin_)*xLen;
	XYZ tAInB = global2local(csB, csA.origin_)*(-xLen);

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
	int spIndexA = it->second;

	it = sphereKeyMap.find(idX2*160000+idY2*400+idZ2);
	if(it == sphereKeyMap.end()){
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

	int indexLen = (int)(len*33.3333);
	int indexAng = (int)(ang*1.111111);
	if(indexAng == 400) indexAng = 0;

	if(indexLen <0 || indexLen >= 500) {
		cout << "invalid indexLen: " << indexLen << endl;
		exit(0);
	}

	if(indexAng < 0 || indexAng >= 400) {
		cout << "invalid index ang: " << indexAng << endl;
		exit(0);
	}

	int index = indexLen*1600 + indexAng*4;

	outIndex[0] = this->biIndex[index];
	outIndex[1] = this->biIndex[index+1];
	outIndex[2] = this->biIndex[index+2];
	outIndex[3] = this->biIndex[index+3];

	outIndex[4] = spIndexA*this->spherePointNum+spIndexB;

	outWeights[0] = this->biWt[index];
	outWeights[1] = this->biWt[index+1];
	outWeights[2] = this->biWt[index+2];
	outWeights[3] = this->biWt[index+3];

	/*
	printf("len interpolation: %6.3f %7.4f %6.3f %6.3f\n", len, indexLen*0.03+0.015, (outIndex[0]/45)*0.3, (outIndex[2]/45)*0.3);
	printf("ang interpolation: %6.2f %7.3f %6.2f %6.2f\n", ang, indexAng*0.8+0.4, (outIndex[0]%45)*8.0, (outIndex[1]%45)*8.0);
	printf("wts: %6.4f %6.4f %6.4f %6.4f\n", outWeights[0], outWeights[1], outWeights[2], outWeights[3]);
	*/
}

void CsMoveTo6DKeySP500::getIndexAndWeight6DInterpolation(const LocalFrame& csA, const LocalFrame& csB, double len, int outIndex[13], double outWeights[13]){


	if(len == 0){
		int index = 0;
		int idA = 0;
		int idB = 0;

		outIndex[0] = this->biIndex[index];
		outIndex[1] = this->biIndex[index+1];
		outIndex[2] = this->biIndex[index+2];
		outIndex[3] = this->biIndex[index+3];

		outIndex[4] = spIndex1[idA]*this->spherePointNum+spIndex1[idB];
		outIndex[5] = spIndex1[idA]*this->spherePointNum+spIndex2[idB];
		outIndex[6] = spIndex1[idA]*this->spherePointNum+spIndex3[idB];
		outIndex[7] = spIndex2[idA]*this->spherePointNum+spIndex1[idB];
		outIndex[8] = spIndex2[idA]*this->spherePointNum+spIndex2[idB];
		outIndex[9] = spIndex2[idA]*this->spherePointNum+spIndex3[idB];
		outIndex[10] = spIndex3[idA]*this->spherePointNum+spIndex1[idB];
		outIndex[11] = spIndex3[idA]*this->spherePointNum+spIndex2[idB];
		outIndex[12] = spIndex3[idA]*this->spherePointNum+spIndex3[idB];


		outWeights[0] = this->biWt[index];
		outWeights[1] = this->biWt[index+1];
		outWeights[2] = this->biWt[index+2];
		outWeights[3] = this->biWt[index+3];
		outWeights[4] = spWt1[idA]*spWt1[idB];
		outWeights[5] = spWt1[idA]*spWt2[idB];
		outWeights[6] = spWt1[idA]*spWt3[idB];
		outWeights[7] = spWt2[idA]*spWt1[idB];
		outWeights[8] = spWt2[idA]*spWt2[idB];
		outWeights[9] = spWt2[idA]*spWt3[idB];
		outWeights[10] = spWt3[idA]*spWt1[idB];
		outWeights[11] = spWt3[idA]*spWt2[idB];
		outWeights[12] = spWt3[idA]*spWt3[idB];

		return;
	}

	double xLen = 1.0/len;

	XYZ tBInA = global2local(csA, csB.origin_)*xLen;
	XYZ tAInB = global2local(csB, csA.origin_)*(-xLen);

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

	it = keyIdToListID.find(idX1*160000+idY1*400+idZ1);
	if(it == keyIdToListID.end()){
		cout << "len: " << len << endl;
		cout << "B in A: " << tBInA.toString() << endl;
		cout << "invalid sphere index: " << idX1 << " " << idY1 << " " << idZ1 << endl;
		exit(1);
	}
	int idA = it->second;

	it = keyIdToListID.find(idX2*160000+idY2*400+idZ2);
	if(it == keyIdToListID.end()){
		cout << "len: " << len << endl;
		cout << "A in B: " << tAInB.toString() << endl;
		cout << "invalid sphere index: " << idX2 << " " << idY2 << " " << idZ2 << endl;
		exit(1);
	}
	int idB = it->second;

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

	int indexLen = (int)(len*33.3333);
	int indexAng = (int)(ang*1.111111);
	if(indexAng == 400) indexAng = 0;

	if(indexLen <0 || indexLen >= 500) {
		cout << "invalid indexLen: " << indexLen << endl;
		exit(0);
	}

	if(indexAng < 0 || indexAng >= 400) {
		cout << "invalid index ang: " << indexAng << endl;
		exit(0);
	}

	int index = indexLen*1600 + indexAng*4;

	outIndex[0] = this->biIndex[index];
	outIndex[1] = this->biIndex[index+1];
	outIndex[2] = this->biIndex[index+2];
	outIndex[3] = this->biIndex[index+3];


	outIndex[4] = spIndex1[idA]*this->spherePointNum+spIndex1[idB];
	outIndex[5] = spIndex1[idA]*this->spherePointNum+spIndex2[idB];
	outIndex[6] = spIndex1[idA]*this->spherePointNum+spIndex3[idB];
	outIndex[7] = spIndex2[idA]*this->spherePointNum+spIndex1[idB];
	outIndex[8] = spIndex2[idA]*this->spherePointNum+spIndex2[idB];
	outIndex[9] = spIndex2[idA]*this->spherePointNum+spIndex3[idB];
	outIndex[10] = spIndex3[idA]*this->spherePointNum+spIndex1[idB];
	outIndex[11] = spIndex3[idA]*this->spherePointNum+spIndex2[idB];
	outIndex[12] = spIndex3[idA]*this->spherePointNum+spIndex3[idB];


	outWeights[0] = this->biWt[index];
	outWeights[1] = this->biWt[index+1];
	outWeights[2] = this->biWt[index+2];
	outWeights[3] = this->biWt[index+3];
	outWeights[4] = spWt1[idA]*spWt1[idB];
	outWeights[5] = spWt1[idA]*spWt2[idB];
	outWeights[6] = spWt1[idA]*spWt3[idB];
	outWeights[7] = spWt2[idA]*spWt1[idB];
	outWeights[8] = spWt2[idA]*spWt2[idB];
	outWeights[9] = spWt2[idA]*spWt3[idB];
	outWeights[10] = spWt3[idA]*spWt1[idB];
	outWeights[11] = spWt3[idA]*spWt2[idB];
	outWeights[12] = spWt3[idA]*spWt3[idB];

	/*
	printf("len interpolation: %6.3f %7.4f %6.3f %6.3f\n", len, indexLen*0.03+0.015, (outIndex[0]/45)*0.3, (outIndex[2]/45)*0.3);
	printf("ang interpolation: %6.2f %7.3f %6.2f %6.2f\n", ang, indexAng*0.8+0.4, (outIndex[0]%45)*8.0, (outIndex[1]%45)*8.0);
	printf("wts: %6.4f %6.4f %6.4f %6.4f\n", outWeights[0], outWeights[1], outWeights[2], outWeights[3]);
	*/
}

int CsMoveTo6DKeySP500::toIndex(const LocalFrame& csA, const LocalFrame& csB, double len){
	if(len == 0){
		return 0;
	}

	double xLen = 1.0/len;

	XYZ tBInA = global2local(csA, csB.origin_)*xLen;
	XYZ tAInB = global2local(csB, csA.origin_)*(-xLen);

	char ss[7];
	ss[6] = '\0';

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

	int spIndexA = sphereKeyMap[idX1*160000+idY1*400+idZ1];
	int spIndexB = sphereKeyMap[idX2*160000+idY2*400+idZ2];

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


	int indexAng = (int)(ang*0.1111111 + 0.5);
	if(indexAng == 40) indexAng = 0;

	int indexLen = (int)(len*3.33333333333+0.5);
	if(indexLen == 50) indexLen = 49;

    return indexLen*10000000 + indexAng*250000 + spIndexA*500 + spIndexB;
}

CsMoveTo6DKeySP500::~CsMoveTo6DKeySP500() {
	// TODO Auto-generated destructor stub
}

Xtb6dEnergy::Xtb6dEnergy(ForceFieldPara* para) {


	string path = NSPdataio::datapath();


	int index, clusterID;
	double ene;
	ifstream file;
	string augc = "AUGC";
    
    cout << "start read file" << endl;

	for(int i=0;i<4;i++){
		for(int j=i;j<4;j++){
			string pairType = augc.substr(i,1) + augc.substr(j,1);
            string fileName = "/public/home/pengx/briqx/xtb/sp500/stat/lowEnergyPoints/"+pairType+".ene";
			file.open(fileName.c_str(), ios::in);

			if(!file.is_open()) {
				cout << "can't open file " << fileName << endl;
			}
            else
                cout << "open file: " << fileName  << endl;
                
			while(file >> index >> ene){
                //cout << index << " " << ene << endl;
				this->nnbKeysEnergy[(i*4+j)][index] = ene;
			}
			file.close();
		}
	}
}

Xtb6dEnergy::~Xtb6dEnergy() {
	// TODO Auto-generated destructor stub
}


}